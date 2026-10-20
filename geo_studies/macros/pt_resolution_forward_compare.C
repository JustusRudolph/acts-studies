#include <TCanvas.h>
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
#include <vector>

#include "../../Utils/BasicFunctions.h"
#include "../../Utils/Acceptance.h"

/*
 * pT resolution from the inverse-pT residual pT_MC * (1/pT_reco - 1/pT_MC) of matched
 * tracks (trackClassification == 1), compared between samples. The truth comes from the
 * track's majority particle in particles_digitized_selected.root, through particle_index
 * (friend file from global/add_particle_index_to_tracks.C); it has to pass the
 * OT_discs_coverage acceptance (0.5 < pT < 10 GeV, 1.5 < |eta| < 2.5) with at least
 * minHits hits.
 * Resolution: sigma of a Gaussian fit over the full residual range (as in O2), per pT
 * and |eta| bin; the RMS is printed alongside to show the tails.
 * Every figure is made separately for primaries (generation 0) and secondaries
 * (generation > 0) of the majority particle.
 */

namespace PtResForward {
const std::vector<int> colours = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<double> ptBins = {0.5, 0.75, 1., 1.5, 2., 3., 5., 10.};
const std::vector<double> etaBins = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
constexpr double resMax = 0.25;
// indexed by isSecondary: generation 0, generation > 0
const std::vector<TString> origins = {"primaries", "secondaries"};

int pdgOf(TString particle) {
  particle.ToLower();
  if (particle.IsDigit()) return particle.Atoi();
  const std::map<TString, int> names = {{"electron", 11}, {"muon", 13}, {"pion", 211},
                                        {"kaon", 321}, {"proton", 2212}};
  return names.count(particle) ? names.at(particle) : -1;
}

std::pair<double, double> gausSigma(TH1D* h) {
  if (h->GetEntries() < 20) return {h->GetRMS(), h->GetRMSError()};
  TF1 gaus("gaus_res", "gaus", -resMax, resMax);
  gaus.SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  h->Fit(&gaus, "QNRL0");
  return {std::abs(gaus.GetParameter(2)), gaus.GetParError(2)};
}

TGraphErrors* resolutionVs(const std::vector<TH1D*>& hists, const std::vector<double>& bins) {
  auto* g = new TGraphErrors();
  for (size_t i_bin = 0; i_bin < hists.size(); i_bin++) {
    if (hists[i_bin]->GetEntries() < 5) continue;
    const auto [sigma, err] = gausSigma(hists[i_bin]);
    g->AddPoint(0.5 * (bins[i_bin] + bins[i_bin + 1]), sigma);
    g->SetPointError(g->GetN() - 1, 0.5 * (bins[i_bin + 1] - bins[i_bin]), err);
  }
  return g;
}

int binOf(const std::vector<double>& bins, double x) {
  const int b = std::upper_bound(bins.begin(), bins.end(), x) - bins.begin() - 1;
  return b >= 0 && b + 1 < (int)bins.size() ? b : -1;
}
}  // namespace PtResForward

/*
 * samples: sample directories in the ACTSO2 output, with seed_<firstSeed> ... seed_<lastSeed>
 * below them, or the sample directory itself if lastSeed < firstSeed. particle: name
 * (electron/muon/pion/kaon/proton) or PDG code, both charges; all particles if empty.
 */
void pt_resolution_forward_compare(const std::vector<TString> samples,
                                   std::vector<TString> labels = {},
                                   const bool onSTBC = false,
                                   const bool isPbPb = false,
                                   const TString particle = "",
                                   const int firstSeed = 0,
                                   const int lastSeed = -1,
                                   const int minHits = 5,
                                   TString tag = "",
                                   const TString trackFile = "tracksummary_ambi-tracks-merged.root") {
  using namespace PtResForward;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Utils::Acceptance acceptance;
  acceptance.standardAcceptance("OT_discs_coverage").addMinHits(minHits);
  const int pdgSelected = particle.IsNull() ? 0 : pdgOf(particle);
  if (pdgSelected < 0) {
    std::cerr << "Error: Unknown particle \"" << particle << "\"." << std::endl;
    return;
  }
  if (pdgSelected > 0) acceptance.addParticleToAccept(pdgSelected);

  if (labels.empty())
    for (const TString& sample : samples) labels.push_back(gSystem->BaseName(sample));
  if (tag.IsNull()) {
    for (const TString& label : labels) tag += (tag.IsNull() ? "" : "_vs_") + label;
    tag.ReplaceAll(" ", "_");
  }
  if (pdgSelected > 0) tag += Form("_pdg%d", pdgSelected);
  if (isPbPb) tag += "_PbPb";
  tag += Form("_minHits%d", minHits);

  std::vector<TString> titles;
  for (const TString& org : origins)
    titles.push_back(org + Form(", 0.5 < p_{T} < 10 GeV, 1.5 < |#eta| < 2.5, #geq %d hits%s", minHits,
                                pdgSelected > 0 ? Form(", |PDG| = %d", pdgSelected) : ""));
  const TString resTitle = "Inverse resolution";

  // [isSecondary][sample]
  std::vector<std::vector<TH1D*>> hRes(origins.size());
  std::vector<std::vector<TGraphErrors*>> gPt(origins.size()), gEta(origins.size());
  for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
    // [isSecondary][bin]
    std::vector<std::vector<TH1D*>> hPt(origins.size()), hEta(origins.size());
    for (const bool isSecondary : {false, true}) {
      const char* org = origins[isSecondary].Data();
      hRes[isSecondary].push_back(new TH1D(Form("hRes_%s_%u", org, i_sample), titles[isSecondary] + ";" + resTitle + ";Tracks", 500, -resMax, resMax));
      for (unsigned i_bin = 0; i_bin + 1 < ptBins.size(); i_bin++)
        hPt[isSecondary].push_back(new TH1D(Form("hResPt_%s_%u_%u", org, i_sample, i_bin), "", 200, -resMax, resMax));
      for (unsigned i_bin = 0; i_bin + 1 < etaBins.size(); i_bin++)
        hEta[isSecondary].push_back(new TH1D(Form("hResEta_%s_%u_%u", org, i_sample, i_bin), "", 200, -resMax, resMax));
    }

    std::vector<TString> dirs;
    for (int seed = firstSeed; seed <= lastSeed; seed++) dirs.push_back(Form("%s/seed_%d", samples[i_sample].Data(), seed));
    // With no seed directories (i.e. a solo run) just use the base directory to get data
    if (dirs.empty()) dirs.push_back(samples[i_sample]);

    // progress every 10% of the directories, if there are more than 10
    const size_t progressStep = dirs.size() > 10 ? (dirs.size() + 9) / 10 : 0;
    for (size_t i_dir = 0; i_dir < dirs.size(); i_dir++) {
      if (progressStep && i_dir % progressStep == 0)
        printf("%s: %zu / %zu directories (%.0f%%)\n", labels[i_sample].Data(), i_dir, dirs.size(),
               100. * i_dir / dirs.size());
      const TString& dir = dirs[i_dir];
      const TString path = Utils::outputPathBase(onSTBC) + dir + "/";
      std::unique_ptr<TFile> fParticles(TFile::Open(path + "particles_digitized_selected.root"));
      std::unique_ptr<TFile> fTracks(TFile::Open(path + trackFile));
      const TString indexFile = path + TString(trackFile).ReplaceAll(".root", "_particle_index.root");
      TTree* tp = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
      TTree* tt = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
      if (tt && !gSystem->AccessPathName(indexFile)) tt->AddFriend("trackindex", indexFile);
      if (!tp || !tt || !tt->GetBranch("particle_index")) {
        std::cerr << "Warning: Missing particles or indexed tracks in " << dir << ", skipping." << std::endl;
        continue;
      }

      // particles, concatenated over events as particle_index counts them
      std::vector<float> pt, eta, phi;
      std::vector<int> pdg, nHits;
      std::vector<unsigned> generation;
      std::vector<float> *vPt = nullptr, *vEta = nullptr, *vPhi = nullptr;
      std::vector<int> *vPdg = nullptr, *vHits = nullptr;
      std::vector<unsigned>* vGen = nullptr;
      tp->SetBranchAddress("pt", &vPt);
      tp->SetBranchAddress("eta", &vEta);
      tp->SetBranchAddress("phi", &vPhi);
      tp->SetBranchAddress("particle_type", &vPdg);
      tp->SetBranchAddress("number_of_hits", &vHits);
      tp->SetBranchAddress("generation", &vGen);
      for (Long64_t i_event = 0; i_event < tp->GetEntries(); i_event++) {
        tp->GetEntry(i_event);
        pt.insert(pt.end(), vPt->begin(), vPt->end());
        eta.insert(eta.end(), vEta->begin(), vEta->end());
        phi.insert(phi.end(), vPhi->begin(), vPhi->end());
        pdg.insert(pdg.end(), vPdg->begin(), vPdg->end());
        nHits.insert(nHits.end(), vHits->begin(), vHits->end());
        generation.insert(generation.end(), vGen->begin(), vGen->end());
      }

      std::vector<int> *index = nullptr, *cls = nullptr;
      std::vector<bool>* hasFit = nullptr;
      std::vector<float> *theta = nullptr, *qop = nullptr;
      tt->SetBranchAddress("particle_index", &index);
      tt->SetBranchAddress("trackClassification", &cls);
      tt->SetBranchAddress("hasFittedParams", &hasFit);
      tt->SetBranchAddress("eTHETA_fit", &theta);
      tt->SetBranchAddress("eQOP_fit", &qop);
      for (Long64_t i_event = 0; i_event < tt->GetEntries(); i_event++) {
        tt->GetEntry(i_event);
        for (size_t i_track = 0; i_track < index->size(); i_track++) {
          const int j = index->at(i_track);
          if (j < 0 || cls->at(i_track) != 1 || !hasFit->at(i_track)) continue;  // 1: matched
          if (!acceptance.isAccepted(pt[j], eta[j], phi[j], std::abs(pdg[j]), nHits[j])) continue;
          const double res = pt[j] * (std::abs(qop->at(i_track)) / std::sin(theta->at(i_track)) - 1. / pt[j]);
          const bool isSecondary = generation[j] > 0;
          hRes[isSecondary][i_sample]->Fill(res);
          if (int b = binOf(ptBins, pt[j]); b >= 0) hPt[isSecondary][b]->Fill(res);
          if (int b = binOf(etaBins, std::abs(eta[j])); b >= 0) hEta[isSecondary][b]->Fill(res);
        }
      }
    }
    for (const bool isSecondary : {false, true}) {
      gPt[isSecondary].push_back(resolutionVs(hPt[isSecondary], ptBins));
      gEta[isSecondary].push_back(resolutionVs(hEta[isSecondary], etaBins));
      const auto [sigma, err] = gausSigma(hRes[isSecondary][i_sample]);
      printf("%s, %s: %.0f tracks, sigma = %.4f +- %.4f (Gaussian), RMS %.4f\n", labels[i_sample].Data(),
             origins[isSecondary].Data(), hRes[isSecondary][i_sample]->GetEntries(), sigma, err, hRes[isSecondary][i_sample]->GetRMS());
    }
  }

  // ---- plots: residuals (unit area), resolution vs pT and vs |eta|
  const TString base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/geo_studies/";
  const TString figureDir = base + "figures/geo/" + tag + "/pt_resolution/";
  gSystem->mkdir(figureDir, true);
  gSystem->mkdir(base + "histos/geo", true);
  auto draw = [&](auto& objs, const TString& name, const TString& title, const TString& axes, double xMin,
                  double xMax, double yMin, double yMax, const char* option) {
    TCanvas c("c_" + name, name, 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetLogy(yMin > 0);
    gPad->DrawFrame(xMin, yMin, xMax, yMax)->SetTitle(title + ";" + axes);
    TLegend leg(0.55, 0.72, 0.88, 0.88);
    leg.SetBorderSize(0);
    for (unsigned i_sample = 0; i_sample < objs.size(); i_sample++) {
      objs[i_sample]->SetLineColor(colours[i_sample % colours.size()]);
      objs[i_sample]->SetMarkerColor(colours[i_sample % colours.size()]);
      objs[i_sample]->SetMarkerStyle(20);
      objs[i_sample]->Draw(option);
      leg.AddEntry(objs[i_sample], labels[i_sample], "l");
    }
    leg.Draw();
    c.SaveAs(figureDir + name + ".pdf");
  };
  for (const bool isSecondary : {false, true}) {
    const TString& org = origins[isSecondary];
    for (TH1D* h : hRes[isSecondary])
      if (h->Integral() > 0) h->Scale(1. / h->Integral());
    draw(hRes[isSecondary], "pt_residuals_" + org, titles[isSecondary], resTitle + ";Fraction of tracks", -resMax, resMax, 1e-5, 1,
         "HIST SAME");
    draw(gPt[isSecondary], "pt_resolution_vs_pt_" + org, titles[isSecondary], "p_{T}^{MC} (GeV);#sigma(p_{T}) / p_{T}",
         ptBins.front(), ptBins.back(), 0, 0.1, "P SAME");
    draw(gEta[isSecondary], "pt_resolution_vs_eta_" + org, titles[isSecondary], "|#eta^{MC}|;#sigma(p_{T}) / p_{T}",
         etaBins.front(), etaBins.back(), 0, 0.1, "P SAME");
  }

  TFile out(base + "histos/geo/pt_resolution_" + tag + ".root", "RECREATE");
  for (const bool isSecondary : {false, true})
    for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
      hRes[isSecondary][i_sample]->Write();
      gPt[isSecondary][i_sample]->Write(Form("gResVsPt_%s_%u", origins[isSecondary].Data(), i_sample));
      gEta[isSecondary][i_sample]->Write(Form("gResVsEta_%s_%u", origins[isSecondary].Data(), i_sample));
    }
  printf("Figures in %s, histograms in %s\n", figureDir.Data(), out.GetName());
}
