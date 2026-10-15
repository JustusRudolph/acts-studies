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
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "../../Utils/BasicFunctions.h"
#include "../../Utils/Acceptance.h"

/*
 * Relative pT resolution, (pT_fit - pT_true) / pT_true, of matched tracks (no duplicates
 * or fakes) in the OT_discs_coverage acceptance (0.5 < pT < 10 GeV, 1.5 < |eta| < 2.5) with
 * at least minHits measurements, compared between samples. Everything comes from the
 * track summary; particles.root is only read to select a particle species.
 * Resolution: sigma of a Gaussian fit to the core (two iterations within +-2 sigma).
 */

namespace PtResForward {
const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<double> ptBins = {0.5, 0.75, 1., 1.5, 2., 3., 5., 10.};
const std::vector<double> etaBins = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
constexpr double resMax = 0.25;
using Barcode = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;

int pdgOf(TString particle) {
  particle.ToLower();
  if (particle.IsDigit()) return particle.Atoi();
  const std::map<TString, int> names = {{"electron", 11}, {"muon", 13}, {"pion", 211},
                                        {"kaon", 321}, {"proton", 2212}};
  return names.count(particle) ? names.at(particle) : -1;
}

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

TGraphErrors* resolutionVs(const std::vector<TH1D*>& hists, const std::vector<double>& bins) {
  auto* g = new TGraphErrors();
  for (size_t b = 0; b < hists.size(); b++) {
    if (hists[b]->GetEntries() < 5) continue;
    const auto [sigma, err] = coreSigma(hists[b]);
    g->AddPoint(0.5 * (bins[b] + bins[b + 1]), sigma);
    g->SetPointError(g->GetN() - 1, 0.5 * (bins[b + 1] - bins[b]), err);
  }
  return g;
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

  std::vector<TH1D*> hRes;
  std::vector<TGraphErrors*> gPt, gEta;
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes.push_back(new TH1D(Form("hRes_%u", s), title + ";(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true};Tracks",
                            500, -resMax, resMax));
    std::vector<TH1D*> hPt, hEta;
    for (size_t b = 0; b + 1 < ptBins.size(); b++)
      hPt.push_back(new TH1D(Form("hResPt_%u_%zu", s, b), "", 100, -resMax, resMax));
    for (size_t b = 0; b + 1 < etaBins.size(); b++)
      hEta.push_back(new TH1D(Form("hResEta_%u_%zu", s, b), "", 100, -resMax, resMax));

    std::vector<TString> dirs;
    for (int seed = firstSeed; seed <= lastSeed; seed++) dirs.push_back(Form("%s/seed_%d", samples[s].Data(), seed));
    if (dirs.empty()) dirs.push_back(samples[s]);

    for (const TString& dir : dirs) {
      const TString path = Utils::outputPathBase(onSTBC) + dir + "/";
      std::unique_ptr<TFile> fTracks(TFile::Open(path + trackFile));
      TTree* t = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
      if (!t) {
        std::cerr << "Warning: No tracks in " << dir << ", skipping." << std::endl;
        continue;
      }

      // |PDG| per particle, only needed to select a species
      std::map<Barcode, int> pdgs;
      if (pdgSelected > 0) {
        std::unique_ptr<TFile> fParticles(TFile::Open(path + "particles.root"));
        TTree* tp = fParticles ? fParticles->Get<TTree>("particles") : nullptr;
        unsigned ev;
        std::vector<unsigned> *v1 = nullptr, *v2 = nullptr, *p = nullptr, *g = nullptr, *sp = nullptr;
        std::vector<int>* type = nullptr;
        tp->SetBranchAddress("event_id", &ev);
        tp->SetBranchAddress("vertex_primary", &v1);
        tp->SetBranchAddress("vertex_secondary", &v2);
        tp->SetBranchAddress("particle", &p);
        tp->SetBranchAddress("generation", &g);
        tp->SetBranchAddress("sub_particle", &sp);
        tp->SetBranchAddress("particle_type", &type);
        for (Long64_t i = 0; i < tp->GetEntries(); i++) {
          tp->GetEntry(i);
          for (size_t k = 0; k < type->size(); k++)
            pdgs[{ev, v1->at(k), v2->at(k), p->at(k), g->at(k), sp->at(k)}] = std::abs(type->at(k));
        }
      }

      unsigned ev;
      std::vector<unsigned> *nMeas = nullptr, *v1 = nullptr, *v2 = nullptr, *p = nullptr, *g = nullptr, *sp = nullptr;
      std::vector<int>* cls = nullptr;
      std::vector<bool>* hasFit = nullptr;
      std::vector<float> *tPt = nullptr, *tEta = nullptr, *tPhi = nullptr, *theta = nullptr, *qop = nullptr;
      t->SetBranchAddress("event_nr", &ev);
      t->SetBranchAddress("nMeasurements", &nMeas);
      t->SetBranchAddress("majorityParticleId_vertex_primary", &v1);
      t->SetBranchAddress("majorityParticleId_vertex_secondary", &v2);
      t->SetBranchAddress("majorityParticleId_particle", &p);
      t->SetBranchAddress("majorityParticleId_generation", &g);
      t->SetBranchAddress("majorityParticleId_sub_particle", &sp);
      t->SetBranchAddress("trackClassification", &cls);
      t->SetBranchAddress("hasFittedParams", &hasFit);
      t->SetBranchAddress("t_pT", &tPt);
      t->SetBranchAddress("t_eta", &tEta);
      t->SetBranchAddress("t_phi", &tPhi);
      t->SetBranchAddress("eTHETA_fit", &theta);
      t->SetBranchAddress("eQOP_fit", &qop);
      for (Long64_t i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        for (size_t k = 0; k < tPt->size(); k++) {
          if (cls->at(k) != 1 || !hasFit->at(k)) continue;  // 1: matched, 2: duplicate, 3: fake
          int pdg = 0;
          if (pdgSelected > 0) {
            auto it = pdgs.find({ev, v1->at(k), v2->at(k), p->at(k), g->at(k), sp->at(k)});
            pdg = it == pdgs.end() ? 0 : it->second;
          }
          if (!acceptance.isAccepted(tPt->at(k), tEta->at(k), tPhi->at(k), pdg, nMeas->at(k))) continue;
          const double res = (std::sin(theta->at(k)) / std::abs(qop->at(k)) - tPt->at(k)) / tPt->at(k);
          hRes.back()->Fill(res);
          const int bPt = std::upper_bound(ptBins.begin(), ptBins.end(), tPt->at(k)) - ptBins.begin() - 1;
          const int bEta = std::upper_bound(etaBins.begin(), etaBins.end(), std::abs(tEta->at(k))) - etaBins.begin() - 1;
          if (bPt >= 0 && bPt < (int)hPt.size()) hPt[bPt]->Fill(res);
          if (bEta >= 0 && bEta < (int)hEta.size()) hEta[bEta]->Fill(res);
        }
      }
    }
    gPt.push_back(resolutionVs(hPt, ptBins));
    gEta.push_back(resolutionVs(hEta, etaBins));
    const auto [sigma, err] = coreSigma(hRes.back());
    printf("%s: %.0f tracks, sigma(pT)/pT = %.4f +- %.4f (core), RMS %.4f\n", labels[s].Data(),
           hRes.back()->GetEntries(), sigma, err, hRes.back()->GetRMS());
  }

  // ---- plots: residuals (unit area), resolution vs pT and vs |eta|
  const TString base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/geo_studies/";
  const TString figureDir = base + "figures/geo/" + tag + "/pt_resolution/";
  gSystem->mkdir(figureDir, true);
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
  draw(hRes, "pt_residuals", "(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true};Fraction of tracks",
       -resMax, resMax, 1e-5, 1, "HIST SAME");
  draw(gPt, "pt_resolution_vs_pt", "p_{T}^{true} (GeV);#sigma(p_{T}) / p_{T}", ptBins.front(),
       ptBins.back(), 0, 0.1, "P SAME");
  draw(gEta, "pt_resolution_vs_eta", "|#eta^{true}|;#sigma(p_{T}) / p_{T}", etaBins.front(),
       etaBins.back(), 0, 0.1, "P SAME");

  TFile out(base + "histos/geo/pt_resolution_" + tag + ".root", "RECREATE");
  for (unsigned s = 0; s < samples.size(); s++) {
    hRes[s]->Write();
    gPt[s]->Write(Form("gResVsPt_%u", s));
    gEta[s]->Write(Form("gResVsEta_%u", s));
  }
  printf("Figures in %s, histograms in %s\n", figureDir.Data(), out.GetName());
}
