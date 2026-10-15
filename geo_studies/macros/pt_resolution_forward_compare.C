#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "../../Utils/BasicFunctions.h"
#include "../../Utils/Acceptance.h"

/*
 * Relative pT resolution, (pT_reco - pT_true) / pT_true, of reconstructible particles
 * above ptMin_GeV in etaMin < (|)eta(|) < etaMax, compared between samples (e.g. the
 * same events simulated with different geometries).
 *
 * Reconstructible: in particles_digitized_selected.root, i.e. primaries with enough
 * measurements (see digiConfig in ACTSO2's ChainConfig.py). Each is matched to its
 * track in the track summary through the barcode of the track's majority particle,
 * keeping only tracks classified as matched (one per particle, no duplicates or fakes).
 * pT_reco comes from the fitted parameters, pT_true and eta from the particle.
 *
 * The resolution is the sigma of a Gaussian fit to the core of the residuals (two
 * iterations within +-2 sigma), overall and in bins of pT and eta.
 */

namespace PtResolution {

const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                 kOrange + 7, kCyan + 2, kViolet + 5};

// pT bin edges (GeV) above ptMin, and the eta bin width
const std::vector<double> ptEdges = {0.5, 0.75, 1., 1.5, 2., 3., 5., 10.};
constexpr double etaBinWidth = 0.05;
// residual histograms: range and bins, coarser for the pT and eta bins with fewer entries
constexpr double resMax = 0.25;
constexpr int nResBins = 500;
constexpr int nResBinsPerBin = 100;

/*
 * PDG code (compared as |pdg|) of a particle name or number, 0 for "" (all particles),
 * -1 if unknown.
 */
int pdgOf(TString particle) {
  particle.ToLower();
  if (particle.IsNull()) return 0;
  if (particle.IsDigit()) return particle.Atoi();
  const std::map<TString, int> names = {
    {"e", 11},     {"electron", 11}, {"mu", 13},  {"muon", 13},   {"pi", 211},
    {"pion", 211}, {"k", 321},       {"kaon", 321}, {"p", 2212}, {"proton", 2212}};
  auto it = names.find(particle);
  return it == names.end() ? -1 : it->second;
}

// sigma (and its error) of a Gaussian fit to the core of h, the RMS if too few entries.
// A likelihood fit, as a bin often holds only O(100) residuals.
std::pair<double, double> coreSigma(TH1D* h) {
  if (h->GetEntries() < 20) return {h->GetRMS(), h->GetRMSError()};
  double mean = h->GetMean(), sigma = h->GetRMS();
  TF1 gaus("gaus_core", "gaus", -resMax, resMax);
  for (int iteration = 0; iteration < 2; iteration++) {
    gaus.SetParameters(h->GetMaximum(), mean, sigma);
    h->Fit(&gaus, "QNRL0", "", mean - 2 * sigma, mean + 2 * sigma);
    mean = gaus.GetParameter(1);
    sigma = std::abs(gaus.GetParameter(2));
  }
  return {sigma, gaus.GetParError(2)};
}

}  // namespace PtResolution

/*
 * samples: sample directories in the ACTSO2 output, each with seed_<firstSeed> ...
 * seed_<lastSeed> below it, or the sample directory itself if lastSeed < firstSeed.
 * labels: legend entries, the sample names if empty. particle: only this particle
 * species (name like "pion"/"kaon"/"proton"/"electron"/"muon" or a PDG code, both
 * charges), all particles if empty. tag: names the figure directory and histogram file,
 * the labels joined by "_vs_" if empty.
 */
void pt_resolution(const std::vector<TString> samples,
                   const bool onSTBC = false,
                   const TString particle = "",
                   TString tag = "",
                   const TString trackFile = "tracksummary_ambi-tracks-merged.root") {
  using namespace PtResolution;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  const int pdgSelected = pdgOf(particle);
  if (pdgSelected < 0) {
    std::cerr << "Error: Unknown particle \"" << particle << "\", give a name (e.g. pion) "
              << "or a PDG code." << std::endl;
    return;
  }
  if (labels.empty())
    for (const TString& sample : samples) labels.push_back(gSystem->BaseName(sample));
  if (labels.size() != samples.size()) {
    std::cerr << "Error: " << labels.size() << " labels for " << samples.size()
              << " samples." << std::endl;
    return;
  }
  if (tag.IsNull()) {
    for (const TString& label : labels) tag += (tag.IsNull() ? "" : "_vs_") + label;
    tag.ReplaceAll(" ", "_");
  }
  if (pdgSelected > 0) tag += Form("_pdg%d", pdgSelected);

  const TString actso2_output_base = Utils::outputPathBase(onSTBC);
  const TString geo_studies_base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/geo_studies/";
  const TString figureDir = geo_studies_base + "figures/geo/" + tag + "/pt_resolution/";
  gSystem->mkdir(figureDir, true);

  // pT bins from ptMin on
  std::vector<double> ptBins = {ptMin_GeV};
  for (double edge : ptEdges)
    if (edge > ptMin_GeV) ptBins.push_back(edge);
  const int nPtBins = ptBins.size() - 1;
  const int nEtaBins = std::lround((etaMax - etaMin) / etaBinWidth);

  const TString selection =
    Form("p_{T} > %.2g GeV, %.2g < %s < %.2g%s", ptMin_GeV, etaMin, absEta ? "|#eta|" : "#eta",
         etaMax, pdgSelected > 0 ? Form(", |PDG| = %d", pdgSelected) : "");

  // per sample: all residuals, and per pT and eta bin
  std::vector<TH1D*> hRes;
  std::vector<std::vector<TH1D*>> hResPt, hResEta;
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes.push_back(new TH1D(Form("hRes_%u", s),
                            selection + ";(p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true};Particles",
                            nResBins, -resMax, resMax));
    hResPt.emplace_back();
    hResEta.emplace_back();
    for (int b = 0; b < nPtBins; b++)
      hResPt[s].push_back(new TH1D(Form("hResPt_%u_%d", s, b), "", nResBinsPerBin, -resMax, resMax));
    for (int b = 0; b < nEtaBins; b++)
      hResEta[s].push_back(new TH1D(Form("hResEta_%u_%d", s, b), "", nResBinsPerBin, -resMax, resMax));
  }

  for (unsigned s = 0; s < samples.size(); s++) {
    std::vector<TString> dirs;
    if (lastSeed < firstSeed)
      dirs.push_back(samples[s]);
    else
      for (int seed = firstSeed; seed <= lastSeed; seed++)
        dirs.push_back(Form("%s/seed_%d", samples[s].Data(), seed));

    unsigned nFilesUsed = 0;
    long nReconstructible = 0, nMatched = 0;
    for (const TString& dir : dirs) {
      TFile* fParticles = TFile::Open(actso2_output_base + dir + "/particles_digitized_selected.root");
      TFile* fTracks = TFile::Open(actso2_output_base + dir + "/" + trackFile);
      TTree* tParticles = fParticles && !fParticles->IsZombie() ? (TTree*)fParticles->Get("particles") : nullptr;
      TTree* tTracks = fTracks && !fTracks->IsZombie() ? (TTree*)fTracks->Get("tracksummary") : nullptr;
      if (!tParticles || !tTracks) {
        std::cerr << "Warning: Missing particles or tracks in " << dir << ", skipping." << std::endl;
        if (fParticles) fParticles->Close();
        if (fTracks) fTracks->Close();
        continue;
      }
      nFilesUsed++;

      // reconstructible particles in the selection, by event and barcode
      unsigned eventId;
      std::vector<unsigned> *vtxPrim = nullptr, *vtxSec = nullptr, *part = nullptr,
                            *gen = nullptr, *subPart = nullptr;
      std::vector<int>* pdg = nullptr;
      std::vector<float> *pt = nullptr, *eta = nullptr;
      tParticles->SetBranchAddress("event_id", &eventId);
      tParticles->SetBranchAddress("vertex_primary", &vtxPrim);
      tParticles->SetBranchAddress("vertex_secondary", &vtxSec);
      tParticles->SetBranchAddress("particle", &part);
      tParticles->SetBranchAddress("generation", &gen);
      tParticles->SetBranchAddress("sub_particle", &subPart);
      tParticles->SetBranchAddress("particle_type", &pdg);
      tParticles->SetBranchAddress("pt", &pt);
      tParticles->SetBranchAddress("eta", &eta);
      for (Long64_t i = 0; i < tParticles->GetEntries(); i++) {
        tParticles->GetEntry(i);
        for (size_t k = 0; k < pt->size(); k++) {
          const double etaSel = absEta ? std::abs(eta->at(k)) : eta->at(k);
          if (pt->at(k) < ptMin_GeV || etaSel < etaMin || etaSel > etaMax) continue;
          if (pdgSelected > 0 && std::abs(pdg->at(k)) != pdgSelected) continue;
          reconstructible[{eventId, vtxPrim->at(k), vtxSec->at(k), part->at(k), gen->at(k),
                           subPart->at(k)}] = {pt->at(k), eta->at(k)};
        }
      }
      nReconstructible += reconstructible.size();

      // their matched tracks
      unsigned eventNr;
      std::vector<unsigned> *tVtxPrim = nullptr, *tVtxSec = nullptr, *tPart = nullptr,
                            *tGen = nullptr, *tSubPart = nullptr;
      std::vector<int>* classification = nullptr;
      std::vector<bool>* hasFit = nullptr;
      std::vector<float> *qop = nullptr, *theta = nullptr;
      tTracks->SetBranchAddress("event_nr", &eventNr);
      tTracks->SetBranchAddress("majorityParticleId_vertex_primary", &tVtxPrim);
      tTracks->SetBranchAddress("majorityParticleId_vertex_secondary", &tVtxSec);
      tTracks->SetBranchAddress("majorityParticleId_particle", &tPart);
      tTracks->SetBranchAddress("majorityParticleId_generation", &tGen);
      tTracks->SetBranchAddress("majorityParticleId_sub_particle", &tSubPart);
      tTracks->SetBranchAddress("trackClassification", &classification);
      tTracks->SetBranchAddress("hasFittedParams", &hasFit);
      tTracks->SetBranchAddress("eQOP_fit", &qop);
      tTracks->SetBranchAddress("eTHETA_fit", &theta);
      for (Long64_t i = 0; i < tTracks->GetEntries(); i++) {
        tTracks->GetEntry(i);
        for (size_t k = 0; k < qop->size(); k++) {
          if (classification->at(k) != 1 || !hasFit->at(k) || qop->at(k) == 0) continue;  // 1: matched
          auto it = reconstructible.find({eventNr, tVtxPrim->at(k), tVtxSec->at(k), tPart->at(k),
                                          tGen->at(k), tSubPart->at(k)});
          if (it == reconstructible.end()) continue;
          const double ptReco = std::sin(theta->at(k)) / std::abs(qop->at(k));
          const double res = (ptReco - it->second.pt) / it->second.pt;
          nMatched++;
          hRes[s]->Fill(res);
          const auto ptBin = std::upper_bound(ptBins.begin(), ptBins.end(), it->second.pt) - ptBins.begin() - 1;
          if (ptBin >= 0 && ptBin < nPtBins) hResPt[s][ptBin]->Fill(res);
          const double etaSel = absEta ? std::abs(it->second.eta) : it->second.eta;
          const int etaBin = std::floor((etaSel - etaMin) / etaBinWidth);
          if (etaBin >= 0 && etaBin < nEtaBins) hResEta[s][etaBin]->Fill(res);
        }
      }
      fParticles->Close();
      fTracks->Close();
    }

    const auto [sigma, sigmaErr] = coreSigma(hRes[s]);
    printf("%s: %u of %zu files, %ld reconstructible particles, %ld with a matched track "
           "(%.1f%%); sigma(pT)/pT = %.4f +- %.4f (core fit), RMS %.4f, %.2f%% outside +-%.2f\n",
           labels[s].Data(), nFilesUsed, dirs.size(), nReconstructible, nMatched,
           nReconstructible ? 100. * nMatched / nReconstructible : 0., sigma, sigmaErr,
           hRes[s]->GetRMS(),
           hRes[s]->GetEntries() ? 100. * (hRes[s]->GetBinContent(0) + hRes[s]->GetBinContent(nResBins + 1)) /
                                     hRes[s]->GetEntries() : 0.,
           resMax);
  }

  // resolution vs pT and vs eta
  std::vector<TGraphErrors*> gPt, gEta;
  for (unsigned s = 0; s < samples.size(); s++) {
    gPt.push_back(new TGraphErrors());
    gEta.push_back(new TGraphErrors());
    for (int b = 0; b < nPtBins; b++) {
      if (hResPt[s][b]->GetEntries() < 5) continue;
      const auto [sigma, sigmaErr] = coreSigma(hResPt[s][b]);
      const int n = gPt[s]->GetN();
      gPt[s]->SetPoint(n, 0.5 * (ptBins[b] + ptBins[b + 1]), sigma);
      gPt[s]->SetPointError(n, 0.5 * (ptBins[b + 1] - ptBins[b]), sigmaErr);
    }
    for (int b = 0; b < nEtaBins; b++) {
      if (hResEta[s][b]->GetEntries() < 5) continue;
      const auto [sigma, sigmaErr] = coreSigma(hResEta[s][b]);
      const int n = gEta[s]->GetN();
      gEta[s]->SetPoint(n, etaMin + (b + 0.5) * etaBinWidth, sigma);
      gEta[s]->SetPointError(n, 0.5 * etaBinWidth, sigmaErr);
    }
  }

  auto style = [](auto* obj, unsigned s) {
    obj->SetLineColor(colors[s % colors.size()]);
    obj->SetMarkerColor(colors[s % colors.size()]);
    obj->SetMarkerStyle(20);
    obj->SetMarkerSize(0.7);
  };
  auto legend = [&](const std::vector<TObject*>& objs, const char* option) {
    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (unsigned s = 0; s < objs.size(); s++) leg->AddEntry(objs[s], labels[s], option);
    leg->Draw();
  };

  // ---- residual distributions, normalised to unit area
  {
    TCanvas* c = new TCanvas("c_res", "pT residuals", 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetLogy();
    double yMax = 0;
    for (TH1D* h : hRes) {
      if (h->Integral() > 0) h->Scale(1. / h->Integral());
      yMax = std::max(yMax, h->GetMaximum());
    }
    for (unsigned s = 0; s < hRes.size(); s++) {
      style(hRes[s], s);
      hRes[s]->GetYaxis()->SetTitle("Fraction of particles");
      hRes[s]->SetMaximum(3 * yMax);
      hRes[s]->SetMinimum(1e-5);
      hRes[s]->Draw(s == 0 ? "HIST" : "HIST SAME");
    }
    legend(std::vector<TObject*>(hRes.begin(), hRes.end()), "l");
    c->SaveAs(figureDir + "pt_residuals.pdf");
  }

  // ---- resolution vs pT and vs eta
  auto drawGraphs = [&](std::vector<TGraphErrors*>& graphs, const TString& name,
                        const TString& xTitle, double xMin, double xMax) {
    TCanvas* c = new TCanvas("c_" + name, name, 800, 600);
    gPad->SetLeftMargin(0.13);
    double yMax = 0;
    for (TGraphErrors* g : graphs)
      for (int i = 0; i < g->GetN(); i++) yMax = std::max(yMax, g->GetY()[i] + g->GetEY()[i]);
    TH1* frame = gPad->DrawFrame(xMin, 0, xMax, 1.3 * (yMax > 0 ? yMax : 0.1));
    frame->SetTitle(selection + ";" + xTitle + ";#sigma(p_{T}) / p_{T}");
    for (unsigned s = 0; s < graphs.size(); s++) {
      style(graphs[s], s);
      graphs[s]->Draw("P SAME");
    }
    legend(std::vector<TObject*>(graphs.begin(), graphs.end()), "lep");
    c->SaveAs(figureDir + name + ".pdf");
  };
  drawGraphs(gPt, "pt_resolution_vs_pt", "p_{T}^{true} (GeV)", ptBins.front(), ptBins.back());
  drawGraphs(gEta, "pt_resolution_vs_eta", absEta ? "|#eta^{true}|" : "#eta^{true}", etaMin, etaMax);

  // ---- histograms, to compare more later
  const TString histFile = geo_studies_base + "histos/geo/pt_resolution_" + tag + ".root";
  gSystem->mkdir(gSystem->DirName(histFile), true);
  TFile out(histFile, "RECREATE");
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes[s]->Write();
    gPt[s]->Write(Form("gResVsPt_%u", s));
    gEta[s]->Write(Form("gResVsEta_%u", s));
  }
  out.Close();
  printf("Figures in %s, histograms in %s\n", figureDir.Data(), histFile.Data());
}
