#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "../../Utils/BasicFunctions.h"
#include "../../Utils/Acceptance.h"

/*
 * Low-pT comparison of tracking efficiency and pT resolution between samples, per particle
 * species (e, pi, p, K, both charges), in linear pT bins of ptStep from ptMin to ptMax.
 * Cuts on the truth particle: absEtaMin < |eta| < absEtaMax, at least minHits hits.
 *
 * Efficiency: denominator are the accepted particles in particles_matched.root (from
 * global/add_matching_to_particles.C), numerator those flagged matched (tracking) or
 * seed_matched (seeding).
 * Resolution: inverse-pT residual pT_MC * (1/pT_reco - 1/pT_MC) of matched tracks
 * (trackClassification == 1), truth from the track's majority particle, looked up in
 * particles_digitized_selected.root by its barcode. Sigma of a Gaussian fit over +-resMax per
 * pT bin (RMS below 20 entries); the pT-integrated residuals are drawn too, to judge the tails.
 * Every efficiency/resolution figure is made separately for primaries (generation 0) and
 * secondaries (generation > 0).
 * Fake rates, vs reconstructed pT with the |eta| cut on the reconstructed eta (fakes have no
 * truth, so no species or origin split):
 *  - seeds: fraction of seeds in iter_<n>/[iter_<n>_]estimatedparameters.root (all iterations,
 *    prevtx excluded) with truthMatched == false, i.e. not all three hits from one particle.
 *  - tracks: fraction of tracks in trackFile with >= minHits measurements that are neither
 *    matched (1) nor duplicate (2). ACTSO2's TrackTruthMatcher leaves tracks whose majority
 *    particle is not in the selected particle collection at 0 (unknown); those count as
 *    unmatched too. The split 0 / 3 (fake) is printed.
 * Output: figures/lowPtComparison/<tag>/ and hists/lowPtComparison/<tag>.root.
 */

namespace LowPtCompare {
const std::vector<int> colours = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<int> pdgCodes = {11, 211, 2212, 321};
const std::vector<TString> names = {"electron", "pion", "proton", "kaon"};
const std::vector<TString> symbols = {"e^{#pm}", "#pi^{#pm}", "p^{#pm}", "K^{#pm}"};
// indexed by isSecondary: generation 0, generation > 0
const std::vector<TString> origins = {"primaries", "secondaries"};

int speciesOf(int pdg) {
  const auto it = std::find(pdgCodes.begin(), pdgCodes.end(), std::abs(pdg));
  return it == pdgCodes.end() ? -1 : it - pdgCodes.begin();
}

std::pair<double, double> gausSigma(TH1D* h, double range) {
  if (h->GetEntries() < 20) return {h->GetRMS(), h->GetRMSError()};
  TF1 gaus("gaus_res", "gaus", -range, range);
  gaus.SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  h->Fit(&gaus, "QNRL0");
  return {std::abs(gaus.GetParameter(2)), gaus.GetParError(2)};
}
}  // namespace LowPtCompare

/*
 * samples: sample directories in the ACTSO2 output, with seed_<firstSeed> ... seed_<lastSeed>
 * below them, or the sample directory itself if lastSeed < firstSeed.
 */
void lowpt_eff_ptres_compare(const std::vector<TString> samples,
                             std::vector<TString> labels = {},
                             const bool onSTBC = false,
                             const int firstSeed = 0,
                             const int lastSeed = -1,
                             const double absEtaMin = 0.,
                             const double absEtaMax = 1.,
                             const int minHits = 5,
                             const double ptMin = 0.,
                             const double ptMax = 0.2,
                             const double ptStep = 0.02,
                             const double resMax = 1.,
                             TString tag = "",
                             const TString trackFile = "tracksummary_ambi-tracks-merged.root") {
  using namespace LowPtCompare;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Utils::Acceptance acceptance;
  acceptance.addPtRange(ptMin, ptMax).addAbsEtaRange(absEtaMin, absEtaMax).addMinHits(minHits);

  if (labels.empty())
    for (const TString& sample : samples) labels.push_back(gSystem->BaseName(sample));
  if (tag.IsNull()) {
    for (const TString& label : labels) tag += (tag.IsNull() ? "" : "_vs_") + label;
    tag.ReplaceAll(" ", "_");
  }
  tag += Form("_eta%.1f-%.1f_minHits%d", absEtaMin, absEtaMax, minHits);

  const int nPtBins = std::lround((ptMax - ptMin) / ptStep);
  const TString selection = Form("%.1f < |#eta| < %.1f, #geq %d hits", absEtaMin, absEtaMax, minHits);

  // [isSecondary][sample][species]
  using PerSpecies = std::vector<std::vector<std::vector<TEfficiency*>>>;
  PerSpecies eff(origins.size(), std::vector<std::vector<TEfficiency*>>(samples.size()));
  PerSpecies effSeed(origins.size(), std::vector<std::vector<TEfficiency*>>(samples.size()));
  std::vector<std::vector<std::vector<TGraphErrors*>>> gRes(origins.size(), std::vector<std::vector<TGraphErrors*>>(samples.size()));
  std::vector<std::vector<std::vector<TH1D*>>> hResAll(origins.size(), std::vector<std::vector<TH1D*>>(samples.size()));
  // [sample], passed = fake
  std::vector<TEfficiency*> fakeSeed, fakeTrack;

  for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
    fakeSeed.push_back(new TEfficiency(Form("fakeSeed_%u", i_sample), "", nPtBins, ptMin, ptMin + nPtBins * ptStep));
    fakeTrack.push_back(new TEfficiency(Form("fakeTrack_%u", i_sample), "", nPtBins, ptMin, ptMin + nPtBins * ptStep));
    // tracks per trackClassification 0..3, within the fake-rate selection
    std::vector<long> nTracksPerClass(4, 0);
    // [isSecondary][species][pT bin]
    std::vector<std::vector<std::vector<TH1D*>>> hRes(origins.size(), std::vector<std::vector<TH1D*>>(names.size()));
    for (const bool isSecondary : {false, true}) {
      const char* org = origins[isSecondary].Data();
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        const char* name = names[i_part].Data();
        eff[isSecondary][i_sample].push_back(
            new TEfficiency(Form("eff_%s_%s_%u", name, org, i_sample), "", nPtBins, ptMin, ptMin + nPtBins * ptStep));
        effSeed[isSecondary][i_sample].push_back(
            new TEfficiency(Form("effSeed_%s_%s_%u", name, org, i_sample), "", nPtBins, ptMin, ptMin + nPtBins * ptStep));
        hResAll[isSecondary][i_sample].push_back(new TH1D(Form("hRes_%s_%s_%u", name, org, i_sample), "", 200, -resMax, resMax));
        for (int i_bin = 0; i_bin < nPtBins; i_bin++)
          hRes[isSecondary][i_part].push_back(new TH1D(Form("hResPt_%s_%s_%u_%d", name, org, i_sample, i_bin), "", 400, -resMax, resMax));
      }
    }

    std::vector<TString> dirs;
    for (int seed = firstSeed; seed <= lastSeed; seed++) dirs.push_back(Form("%s/seed_%d", samples[i_sample].Data(), seed));
    // With no seed directories (i.e. a solo run) just use the base directory to get data
    if (dirs.empty()) dirs.push_back(samples[i_sample]);

    // progress every 10% of the directories, if there are more than 10
    const size_t progressStep = dirs.size() > 10 ? (dirs.size() + 9) / 10 : 0;
    for (size_t i_dir = 0; i_dir < dirs.size(); i_dir++) {
      if (progressStep && i_dir % progressStep == 0)
        printf("%s: %zu / %zu directories (%.0f%%)\n", labels[i_sample].Data(), i_dir, dirs.size(), 100. * i_dir / dirs.size());
      const TString path = Utils::outputPathBase(onSTBC) + dirs[i_dir] + "/";

      // ---- efficiency
      {
        std::unique_ptr<TFile> fMatched(TFile::Open(path + "particles_matched.root"));
        TTree* tp = fMatched && !fMatched->IsZombie() ? fMatched->Get<TTree>("particles") : nullptr;
        if (!tp) {
          std::cerr << "Warning: Cannot read " << path << "particles_matched.root, no efficiency from it." << std::endl;
        } else {
          std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr;
          std::vector<int> *pdg = nullptr, *nHits = nullptr;
          std::vector<unsigned>* generation = nullptr;
          std::vector<bool>* matched = nullptr;
          std::vector<bool>* seedMatched = nullptr;
          tp->SetBranchAddress("pt", &pt);
          tp->SetBranchAddress("eta", &eta);
          tp->SetBranchAddress("phi", &phi);
          tp->SetBranchAddress("particle_type", &pdg);
          tp->SetBranchAddress("number_of_hits", &nHits);
          tp->SetBranchAddress("generation", &generation);
          tp->SetBranchAddress("matched", &matched);
          if (tp->GetBranch("seed_matched"))
            tp->SetBranchAddress("seed_matched", &seedMatched);
          else
            std::cerr << "Warning: No seed_matched in " << path << "particles_matched.root, no seeding efficiency from it." << std::endl;
          for (Long64_t i_event = 0; i_event < tp->GetEntries(); i_event++) {
            tp->GetEntry(i_event);
            for (size_t i_particle = 0; i_particle < pt->size(); i_particle++) {
              const int i_part = speciesOf(pdg->at(i_particle));
              if (i_part < 0) continue;
              if (!acceptance.isAccepted(pt->at(i_particle), eta->at(i_particle), phi->at(i_particle),
                                         std::abs(pdg->at(i_particle)), nHits->at(i_particle)))
                continue;
              const bool isSecondary = generation->at(i_particle) > 0;
              eff[isSecondary][i_sample][i_part]->Fill(matched->at(i_particle), pt->at(i_particle));
              if (seedMatched) effSeed[isSecondary][i_sample][i_part]->Fill(seedMatched->at(i_particle), pt->at(i_particle));
            }
          }
        }
      }

      // ---- seed fake rate, over all iterations
      for (int i_iter = 0;; i_iter++) {
        const TString iterDir = path + Form("iter_%d/", i_iter);
        if (gSystem->AccessPathName(iterDir)) {
          if (i_iter == 0) std::cerr << "Warning: No " << iterDir << ", no seed fake rate from it." << std::endl;
          break;
        }
        // iteration 0 writes estimatedparameters.root, later ones iter_<n>_estimatedparameters.root
        TString seedFile = iterDir + Form("iter_%d_estimatedparameters.root", i_iter);
        if (gSystem->AccessPathName(seedFile)) seedFile = iterDir + "estimatedparameters.root";
        std::unique_ptr<TFile> fSeeds(TFile::Open(seedFile));
        TTree* ts = fSeeds && !fSeeds->IsZombie() ? fSeeds->Get<TTree>("estimatedparams") : nullptr;
        if (!ts) {
          std::cerr << "Warning: Cannot read " << seedFile << std::endl;
          continue;
        }
        float seedPt, seedEta;
        bool seedTruthMatched;
        ts->SetBranchStatus("*", false);
        for (const char* name : {"pt", "eta", "truthMatched"}) ts->SetBranchStatus(name, true);
        ts->SetBranchAddress("pt", &seedPt);
        ts->SetBranchAddress("eta", &seedEta);
        ts->SetBranchAddress("truthMatched", &seedTruthMatched);
        for (Long64_t i_seed = 0; i_seed < ts->GetEntries(); i_seed++) {
          ts->GetEntry(i_seed);
          if (std::abs(seedEta) < absEtaMin || std::abs(seedEta) >= absEtaMax) continue;
          fakeSeed[i_sample]->Fill(!seedTruthMatched, seedPt);
        }
      }

      // ---- track fake rate and resolution
      std::unique_ptr<TFile> fTracks(TFile::Open(path + trackFile));
      TTree* tt = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
      if (!tt) {
        std::cerr << "Warning: Cannot read " << path << trackFile << ", no track fake rate or resolution from it." << std::endl;
        continue;
      }

      // particles, concatenated over events, and their index by barcode
      // (event, vertex_primary, vertex_secondary, particle, generation, sub_particle)
      using Barcode = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;
      std::map<Barcode, int> indexOf;
      std::vector<float> pt, eta, phi;
      std::vector<int> pdg, nHits;
      std::vector<unsigned> generation;
      std::unique_ptr<TFile> fParticles(TFile::Open(path + "particles_digitized_selected.root"));
      TTree* tp = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
      if (!tp) {
        std::cerr << "Warning: Cannot read " << path << "particles_digitized_selected.root, no resolution from it." << std::endl;
      } else {
        unsigned eventId;
        std::vector<float> *vPt = nullptr, *vEta = nullptr, *vPhi = nullptr;
        std::vector<int> *vPdg = nullptr, *vHits = nullptr;
        std::vector<unsigned> *vGen = nullptr, *vVtx1 = nullptr, *vVtx2 = nullptr, *vPart = nullptr, *vSub = nullptr;
        tp->SetBranchAddress("event_id", &eventId);
        tp->SetBranchAddress("pt", &vPt);
        tp->SetBranchAddress("eta", &vEta);
        tp->SetBranchAddress("phi", &vPhi);
        tp->SetBranchAddress("particle_type", &vPdg);
        tp->SetBranchAddress("number_of_hits", &vHits);
        tp->SetBranchAddress("generation", &vGen);
        tp->SetBranchAddress("vertex_primary", &vVtx1);
        tp->SetBranchAddress("vertex_secondary", &vVtx2);
        tp->SetBranchAddress("particle", &vPart);
        tp->SetBranchAddress("sub_particle", &vSub);
        for (Long64_t i_event = 0; i_event < tp->GetEntries(); i_event++) {
          tp->GetEntry(i_event);
          for (size_t i_particle = 0; i_particle < vPt->size(); i_particle++)
            indexOf[{eventId, vVtx1->at(i_particle), vVtx2->at(i_particle), vPart->at(i_particle), vGen->at(i_particle),
                     vSub->at(i_particle)}] = pt.size() + i_particle;
          pt.insert(pt.end(), vPt->begin(), vPt->end());
          eta.insert(eta.end(), vEta->begin(), vEta->end());
          phi.insert(phi.end(), vPhi->begin(), vPhi->end());
          pdg.insert(pdg.end(), vPdg->begin(), vPdg->end());
          nHits.insert(nHits.end(), vHits->begin(), vHits->end());
          generation.insert(generation.end(), vGen->begin(), vGen->end());
        }
      }

      unsigned eventNr;
      std::vector<int>* cls = nullptr;
      std::vector<unsigned>* nMeas = nullptr;
      std::vector<bool>* hasFit = nullptr;
      std::vector<float> *theta = nullptr, *qop = nullptr;
      std::vector<unsigned> *mVtx1 = nullptr, *mVtx2 = nullptr, *mPart = nullptr, *mGen = nullptr, *mSub = nullptr;
      tt->SetBranchStatus("*", false);
      for (const char* name : {"event_nr", "trackClassification", "nMeasurements", "hasFittedParams", "eTHETA_fit",
                               "eQOP_fit", "majorityParticleId_vertex_primary", "majorityParticleId_vertex_secondary",
                               "majorityParticleId_particle", "majorityParticleId_generation",
                               "majorityParticleId_sub_particle"})
        tt->SetBranchStatus(name, true);
      tt->SetBranchAddress("event_nr", &eventNr);
      tt->SetBranchAddress("trackClassification", &cls);
      tt->SetBranchAddress("nMeasurements", &nMeas);
      tt->SetBranchAddress("hasFittedParams", &hasFit);
      tt->SetBranchAddress("eTHETA_fit", &theta);
      tt->SetBranchAddress("eQOP_fit", &qop);
      tt->SetBranchAddress("majorityParticleId_vertex_primary", &mVtx1);
      tt->SetBranchAddress("majorityParticleId_vertex_secondary", &mVtx2);
      tt->SetBranchAddress("majorityParticleId_particle", &mPart);
      tt->SetBranchAddress("majorityParticleId_generation", &mGen);
      tt->SetBranchAddress("majorityParticleId_sub_particle", &mSub);
      for (Long64_t i_event = 0; i_event < tt->GetEntries(); i_event++) {
        tt->GetEntry(i_event);
        for (size_t i_track = 0; i_track < cls->size(); i_track++) {
          if (!hasFit->at(i_track)) continue;
          const double recoPt = std::abs(std::sin(theta->at(i_track)) / qop->at(i_track));
          const double recoEta = -std::log(std::tan(0.5 * theta->at(i_track)));

          // fake rate: neither matched (1) nor duplicate (2)
          if (std::abs(recoEta) >= absEtaMin && std::abs(recoEta) < absEtaMax && nMeas->at(i_track) >= unsigned(minHits)) {
            const int trackClass = cls->at(i_track);
            fakeTrack[i_sample]->Fill(trackClass != 1 && trackClass != 2, recoPt);
            if (recoPt >= ptMin && recoPt < ptMax && trackClass >= 0 && trackClass < 4) nTracksPerClass[trackClass]++;
          }

          // resolution of matched tracks
          if (cls->at(i_track) != 1) continue;
          const auto it = indexOf.find({eventNr, mVtx1->at(i_track), mVtx2->at(i_track), mPart->at(i_track),
                                        mGen->at(i_track), mSub->at(i_track)});
          if (it == indexOf.end()) continue;
          const int j = it->second;
          const int i_part = speciesOf(pdg[j]);
          if (i_part < 0 || !acceptance.isAccepted(pt[j], eta[j], phi[j], std::abs(pdg[j]), nHits[j])) continue;
          const int i_bin = std::floor((pt[j] - ptMin) / ptStep);
          if (i_bin < 0 || i_bin >= nPtBins) continue;
          const double res = pt[j] / recoPt - 1.;  // pT_MC * (1/pT_reco - 1/pT_MC)
          const bool isSecondary = generation[j] > 0;
          hRes[isSecondary][i_part][i_bin]->Fill(res);
          hResAll[isSecondary][i_sample][i_part]->Fill(res);
        }
      }
    }

    {
      const TH1* total = fakeSeed[i_sample]->GetTotalHistogram();
      const TH1* passed = fakeSeed[i_sample]->GetPassedHistogram();
      printf("%s, %.2g < p_T^reco < %.2g GeV: seed fake rate %.3f (%.0f seeds)\n", labels[i_sample].Data(), ptMin, ptMax,
             total->Integral() ? passed->Integral() / total->Integral() : 0., total->Integral());
      total = fakeTrack[i_sample]->GetTotalHistogram();
      passed = fakeTrack[i_sample]->GetPassedHistogram();
      printf("%s, %.2g < p_T^reco < %.2g GeV: track unmatched rate %.3f (%.0f tracks; unknown %ld, matched %ld, duplicate %ld, fake %ld)\n",
             labels[i_sample].Data(), ptMin, ptMax, total->Integral() ? passed->Integral() / total->Integral() : 0.,
             total->Integral(), nTracksPerClass[0], nTracksPerClass[1], nTracksPerClass[2], nTracksPerClass[3]);
    }

    // resolution vs pT per species, and a summary line per species
    for (const bool isSecondary : {false, true}) {
      printf("%s, %s, %.2g < p_T < %.2g GeV:\n", labels[i_sample].Data(), origins[isSecondary].Data(), ptMin, ptMax);
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        auto* g = new TGraphErrors();
        for (int i_bin = 0; i_bin < nPtBins; i_bin++) {
          TH1D* h = hRes[isSecondary][i_part][i_bin];
          if (h->GetEntries() < 5) continue;
          const auto [sigma, err] = gausSigma(h, resMax);
          g->AddPoint(ptMin + (i_bin + 0.5) * ptStep, sigma);
          g->SetPointError(g->GetN() - 1, 0.5 * ptStep, err);
        }
        gRes[isSecondary][i_sample].push_back(g);
        for (TH1D* h : hRes[isSecondary][i_part]) delete h;

        const TH1* total = eff[isSecondary][i_sample][i_part]->GetTotalHistogram();
        const TH1* passed = eff[isSecondary][i_sample][i_part]->GetPassedHistogram();
        const TH1* totalSeed = effSeed[isSecondary][i_sample][i_part]->GetTotalHistogram();
        const TH1* passedSeed = effSeed[isSecondary][i_sample][i_part]->GetPassedHistogram();
        printf("  %-8s seeding efficiency %.3f, tracking efficiency %.3f (%.0f particles), %.0f matched tracks, residual RMS %.4f\n",
               names[i_part].Data(), totalSeed->GetEntries() ? passedSeed->GetEntries() / totalSeed->GetEntries() : 0.,
               total->GetEntries() ? passed->GetEntries() / total->GetEntries() : 0., total->GetEntries(),
               hResAll[isSecondary][i_sample][i_part]->GetEntries(), hResAll[isSecondary][i_sample][i_part]->GetRMS());
      }
    }
  }

  // ---- plots: one pad per species, samples overlaid
  const TString base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/pythia_studies/";
  const TString figureDir = base + "figures/lowPtComparison/" + tag + "/";
  gSystem->mkdir(figureDir, true);
  gSystem->mkdir(base + "hists/lowPtComparison", true);

  // objs[sample][species]; yMax <= 0: from the data
  auto draw = [&](auto& objs, const TString& name, const TString& title, const TString& axes, double xMin, double xMax,
                  double yMin, double yMax, const char* option, bool logy, double legY) {
    TCanvas c("c_" + name, name, 1200, 900);
    c.Divide(2, 2);
    std::vector<TLegend*> legends;
    for (unsigned i_part = 0; i_part < names.size(); i_part++) {
      c.cd(i_part + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetGridy();
      gPad->SetLogy(logy);
      double padYMax = yMax;
      if (padYMax <= 0) {
        for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
          auto* g = dynamic_cast<TGraphErrors*>(objs[i_sample][i_part]);
          for (int i_point = 0; g && i_point < g->GetN(); i_point++)
            padYMax = std::max(padYMax, g->GetPointY(i_point) + g->GetErrorY(i_point));
        }
        padYMax = padYMax > 0 ? 1.25 * padYMax : 0.1;
      }
      gPad->DrawFrame(xMin, yMin, xMax, padYMax)->SetTitle(symbols[i_part] + " " + title + ";" + axes);
      legends.push_back(new TLegend(0.55, legY - 0.06 * samples.size(), 0.88, legY));
      legends.back()->SetBorderSize(0);
      legends.back()->SetFillStyle(0);
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        auto* obj = objs[i_sample][i_part];
        obj->SetLineColor(colours[i_sample % colours.size()]);
        obj->SetMarkerColor(colours[i_sample % colours.size()]);
        obj->SetMarkerStyle(20);
        obj->SetMarkerSize(0.8);
        obj->Draw(option);
        legends.back()->AddEntry(obj, labels[i_sample], "lp");
      }
      legends.back()->Draw();
    }
    c.SaveAs(figureDir + name + ".pdf");
    for (TLegend* leg : legends) delete leg;
  };

  for (const bool isSecondary : {false, true}) {
    const TString& org = origins[isSecondary];
    const TString title = org + ", " + selection;
    draw(eff[isSecondary], "efficiency_vs_pT_" + org, title, "p_{T} (GeV/c);Efficiency", ptMin, ptMax, 0, 1.05, "P SAME",
         false, 0.35);
    draw(effSeed[isSecondary], "seeding_efficiency_vs_pT_" + org, title, "p_{T} (GeV/c);Seeding efficiency", ptMin, ptMax, 0,
         1.05, "P SAME", false, 0.35);
    draw(gRes[isSecondary], "pt_resolution_vs_pT_" + org, title, "p_{T}^{MC} (GeV/c);#sigma(p_{T}) / p_{T}", ptMin, ptMax, 0,
         -1, "P SAME", false, 0.88);
    for (auto& perSample : hResAll[isSecondary])
      for (TH1D* h : perSample)
        if (h->Integral() > 0) h->Scale(1. / h->Integral());
    draw(hResAll[isSecondary], "pt_residuals_" + org, title + Form(", %.2g < p_{T} < %.2g GeV/c", ptMin, ptMax),
         "p_{T}^{MC} (1/p_{T}^{reco} - 1/p_{T}^{MC});Fraction of tracks", -resMax, resMax, 1e-4, 1, "HIST SAME", true, 0.88);
  }

  // fake rates: one pad, samples overlaid
  auto drawFake = [&](std::vector<TEfficiency*>& objs, const TString& name, const TString& title, const TString& axes) {
    TCanvas c("c_" + name, name, 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetGridy();
    gPad->DrawFrame(ptMin, 0, ptMax, 1.05)->SetTitle(title + ";" + axes);
    TLegend leg(0.55, 0.88 - 0.06 * samples.size(), 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
      objs[i_sample]->SetLineColor(colours[i_sample % colours.size()]);
      objs[i_sample]->SetMarkerColor(colours[i_sample % colours.size()]);
      objs[i_sample]->SetMarkerStyle(20);
      objs[i_sample]->SetMarkerSize(0.8);
      objs[i_sample]->Draw("P SAME");
      leg.AddEntry(objs[i_sample], labels[i_sample], "lp");
    }
    leg.Draw();
    c.SaveAs(figureDir + name + ".pdf");
  };
  const TString recoSelection = Form("%.1f < |#eta^{reco}| < %.1f", absEtaMin, absEtaMax);
  drawFake(fakeSeed, "seed_fake_rate_vs_pT", "Seeds (all iterations), " + recoSelection,
           "p_{T}^{reco} (GeV/c);Fake seed rate");
  drawFake(fakeTrack, "track_fake_rate_vs_pT", "Tracks, " + recoSelection + Form(", #geq %d measurements", minHits),
           "p_{T}^{reco} (GeV/c);Unmatched track rate");

  TFile out(base + "hists/lowPtComparison/" + tag + ".root", "RECREATE");
  for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
    fakeSeed[i_sample]->Write();
    fakeTrack[i_sample]->Write();
  }
  for (const bool isSecondary : {false, true})
    for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++)
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        eff[isSecondary][i_sample][i_part]->Write();
        effSeed[isSecondary][i_sample][i_part]->Write();
        hResAll[isSecondary][i_sample][i_part]->Write();
        gRes[isSecondary][i_sample][i_part]->Write(
            Form("gRes_%s_%s_%u", names[i_part].Data(), origins[isSecondary].Data(), i_sample));
      }
  printf("Figures in %s, histograms in %s\n", figureDir.Data(), out.GetName());
}
