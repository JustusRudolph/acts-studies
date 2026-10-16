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
 * from add_particle_index_to_tracks.C; it has to pass the OT_discs_coverage acceptance
 * (0.5 < pT < 10 GeV, 1.5 < |eta| < 2.5) with at least minHits hits.
 * Resolution: sigma of a Gaussian fit over the full residual range (as in O2), per pT
 * and |eta| bin; the RMS is printed alongside to show the tails.
 */

namespace PtResForward {
const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<double> ptBins = {0.5, 0.75, 1., 1.5, 2., 3., 5., 10.};
const std::vector<double> etaBins = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
constexpr double resMax = 0.25;

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
  for (size_t b = 0; b < hists.size(); b++) {
    if (hists[b]->GetEntries() < 5) continue;
    const auto [sigma, err] = gausSigma(hists[b]);
    g->AddPoint(0.5 * (bins[b] + bins[b + 1]), sigma);
    g->SetPointError(g->GetN() - 1, 0.5 * (bins[b + 1] - bins[b]), err);
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
  const TString title = Form("0.5 < p_{T} < 10 GeV, 1.5 < |#eta| < 2.5, #geq %d hits%s", minHits,
                             pdgSelected > 0 ? Form(", |PDG| = %d", pdgSelected) : "");
  const TString resTitle = "p_{T}^{MC} (1/p_{T}^{reco} - 1/p_{T}^{MC})";

  std::vector<TH1D*> hRes;
  std::vector<TGraphErrors*> gPt, gEta;
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes.push_back(new TH1D(Form("hRes_%u", s), title + ";" + resTitle + ";Tracks", 500, -resMax, resMax));
    std::vector<TH1D*> hPt, hEta;
    for (size_t b = 0; b + 1 < ptBins.size(); b++)
      hPt.push_back(new TH1D(Form("hResPt_%u_%zu", s, b), "", 200, -resMax, resMax));
    for (size_t b = 0; b + 1 < etaBins.size(); b++)
      hEta.push_back(new TH1D(Form("hResEta_%u_%zu", s, b), "", 200, -resMax, resMax));

    std::vector<TString> dirs;
    for (int seed = firstSeed; seed <= lastSeed; seed++) dirs.push_back(Form("%s/seed_%d", samples[s].Data(), seed));
    if (dirs.empty()) dirs.push_back(samples[s]);

    for (const TString& dir : dirs) {
      const TString path = Utils::outputPathBase(onSTBC) + dir + "/";
      std::unique_ptr<TFile> fParticles(TFile::Open(path + "particles_digitized_selected.root"));
      std::unique_ptr<TFile> fTracks(TFile::Open(path + trackFile));
      TTree* tp = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
      TTree* tt = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
      if (!tp || !tt || !tt->GetBranch("particle_index")) {
        std::cerr << "Warning: Missing particles or indexed tracks in " << dir << ", skipping." << std::endl;
        continue;
      }

      // particles, concatenated over events as particle_index counts them
      std::vector<float> pt, eta, phi;
      std::vector<int> pdg, nHits;
      std::vector<float> *vPt = nullptr, *vEta = nullptr, *vPhi = nullptr;
      std::vector<int> *vPdg = nullptr, *vHits = nullptr;
      tp->SetBranchAddress("pt", &vPt);
      tp->SetBranchAddress("eta", &vEta);
      tp->SetBranchAddress("phi", &vPhi);
      tp->SetBranchAddress("particle_type", &vPdg);
      tp->SetBranchAddress("number_of_hits", &vHits);
      for (Long64_t i = 0; i < tp->GetEntries(); i++) {
        tp->GetEntry(i);
        pt.insert(pt.end(), vPt->begin(), vPt->end());
        eta.insert(eta.end(), vEta->begin(), vEta->end());
        phi.insert(phi.end(), vPhi->begin(), vPhi->end());
        pdg.insert(pdg.end(), vPdg->begin(), vPdg->end());
        nHits.insert(nHits.end(), vHits->begin(), vHits->end());
      }

      std::vector<int> *index = nullptr, *cls = nullptr;
      std::vector<bool>* hasFit = nullptr;
      std::vector<float> *theta = nullptr, *qop = nullptr;
      tt->SetBranchAddress("particle_index", &index);
      tt->SetBranchAddress("trackClassification", &cls);
      tt->SetBranchAddress("hasFittedParams", &hasFit);
      tt->SetBranchAddress("eTHETA_fit", &theta);
      tt->SetBranchAddress("eQOP_fit", &qop);
      for (Long64_t i = 0; i < tt->GetEntries(); i++) {
        tt->GetEntry(i);
        for (size_t k = 0; k < index->size(); k++) {
          const int j = index->at(k);
          if (j < 0 || cls->at(k) != 1 || !hasFit->at(k)) continue;  // 1: matched
          if (!acceptance.isAccepted(pt[j], eta[j], phi[j], std::abs(pdg[j]), nHits[j])) continue;
          const double res = pt[j] * (std::abs(qop->at(k)) / std::sin(theta->at(k)) - 1. / pt[j]);
          hRes.back()->Fill(res);
          if (int b = binOf(ptBins, pt[j]); b >= 0) hPt[b]->Fill(res);
          if (int b = binOf(etaBins, std::abs(eta[j])); b >= 0) hEta[b]->Fill(res);
        }
      }
    }
    gPt.push_back(resolutionVs(hPt, ptBins));
    gEta.push_back(resolutionVs(hEta, etaBins));
    const auto [sigma, err] = gausSigma(hRes.back());
    printf("%s: %.0f tracks, sigma = %.4f +- %.4f (Gaussian), RMS %.4f\n", labels[s].Data(),
           hRes.back()->GetEntries(), sigma, err, hRes.back()->GetRMS());
  }

  // ---- plots: residuals (unit area), resolution vs pT and vs |eta|
  const TString base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/geo_studies/";
  const TString figureDir = base + "figures/geo/" + tag + "/pt_resolution/";
  gSystem->mkdir(figureDir, true);
  gSystem->mkdir(base + "histos/geo", true);
  auto draw = [&](auto& objs, const TString& name, const TString& axes, double xMin, double xMax,
                  double yMin, double yMax, const char* option) {
    TCanvas c("c_" + name, name, 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetLogy(yMin > 0);
    gPad->DrawFrame(xMin, yMin, xMax, yMax)->SetTitle(title + ";" + axes);
    TLegend leg(0.55, 0.72, 0.88, 0.88);
    leg.SetBorderSize(0);
    for (unsigned s = 0; s < objs.size(); s++) {
      objs[s]->SetLineColor(colors[s % colors.size()]);
      objs[s]->SetMarkerColor(colors[s % colors.size()]);
      objs[s]->SetMarkerStyle(20);
      objs[s]->Draw(option);
      leg.AddEntry(objs[s], labels[s], "l");
    }
    leg.Draw();
    c.SaveAs(figureDir + name + ".pdf");
  };
  for (TH1D* h : hRes)
    if (h->Integral() > 0) h->Scale(1. / h->Integral());
  draw(hRes, "pt_residuals", resTitle + ";Fraction of tracks", -resMax, resMax, 1e-5, 1, "HIST SAME");
  draw(gPt, "pt_resolution_vs_pt", "p_{T}^{MC} (GeV);#sigma(p_{T}) / p_{T}", ptBins.front(),
       ptBins.back(), 0, 0.1, "P SAME");
  draw(gEta, "pt_resolution_vs_eta", "|#eta^{MC}|;#sigma(p_{T}) / p_{T}", etaBins.front(),
       etaBins.back(), 0, 0.1, "P SAME");

  TFile out(base + "histos/geo/pt_resolution_" + tag + ".root", "RECREATE");
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes[s]->Write();
    gPt[s]->Write(Form("gResVsPt_%u", s));
    gEta[s]->Write(Form("gResVsEta_%u", s));
  }
  printf("Figures in %s, histograms in %s\n", figureDir.Data(), out.GetName());
}
