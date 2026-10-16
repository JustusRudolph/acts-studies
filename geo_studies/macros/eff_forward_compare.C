#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../Utils/BasicFunctions.h"
#include "../../Utils/Acceptance.h"

/*
 * Tracking efficiency per particle species (e, pi, p, K, both charges) in the forward
 * region 1.5 < |eta| < 2.5, compared between samples: vs pT, and vs |eta| for
 * pT > ptCutForEta. Denominator: primaries in particles_matched.root (from
 * global/add_matching_to_particles.C) with at least minHits hits (7, as
 * alice3_reproduce.C). Numerator: those flagged matched.
 */

namespace EffForward {
const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<int> pdgCodes = {11, 211, 2212, 321};
const std::vector<TString> names = {"electron", "pion", "proton", "kaon"};
const std::vector<TString> symbols = {"e^{#pm}", "#pi^{#pm}", "p^{#pm}", "K^{#pm}"};
constexpr double etaMin = 1.5, etaMax = 2.5;
}  // namespace EffForward

/*
 * samples: sample directories in the ACTSO2 output, with seed_<firstSeed> ... seed_<lastSeed>
 * below them, or the sample directory itself if lastSeed < firstSeed.
 */
void eff_forward_compare(const std::vector<TString> samples,
                         std::vector<TString> labels = {},
                         const bool onSTBC = false,
                         const int firstSeed = 0,
                         const int lastSeed = -1,
                         const int minHits = 7,
                         const double ptCutForEta = 0.5,
                         TString tag = "") {
  using namespace EffForward;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Utils::Acceptance acceptance;
  acceptance.addAbsEtaRange(etaMin, etaMax).addMinHits(minHits);

  if (labels.empty())
    for (const TString& sample : samples) labels.push_back(gSystem->BaseName(sample));
  if (tag.IsNull()) {
    for (const TString& label : labels) tag += (tag.IsNull() ? "" : "_vs_") + label;
    tag.ReplaceAll(" ", "_");
  }

  // log pT bins from 10 MeV to 10 GeV, as in alice3_reproduce.C
  const int nPtBins = 20;
  std::vector<double> ptBins(nPtBins + 1);
  for (int i_bin = 0; i_bin <= nPtBins; i_bin++) ptBins[i_bin] = std::pow(10., -2. + 3. * i_bin / nPtBins);

  // [sample][species]
  std::vector<std::vector<TEfficiency*>> effPt(samples.size()), effEta(samples.size());
  for (unsigned s = 0; s < samples.size(); s++) {
    for (unsigned i_part = 0; i_part < names.size(); i_part++) {
      effPt[s].push_back(new TEfficiency(Form("effPt_%s_%u", names[i_part].Data(), s), "",
                                         nPtBins, ptBins.data()));
      effEta[s].push_back(new TEfficiency(Form("effEta_%s_%u", names[i_part].Data(), s), "",
                                          20, etaMin, etaMax));
    }

    std::vector<TString> dirs;
    for (int seed = firstSeed; seed <= lastSeed; seed++) dirs.push_back(Form("%s/seed_%d", samples[s].Data(), seed));
    // With no seed directories (i.e. a solo run) just use the base directory to get data
    if (dirs.empty()) dirs.push_back(samples[s]);

    // progress every 10% of the directories, if there are more than 10
    const size_t progressStep = dirs.size() > 10 ? (dirs.size() + 9) / 10 : 0;
    for (size_t i_dir = 0; i_dir < dirs.size(); i_dir++) {
      if (progressStep && i_dir % progressStep == 0)
        printf("%s: %zu / %zu directories (%.0f%%)\n", labels[s].Data(), i_dir, dirs.size(),
               100. * i_dir / dirs.size());
      const TString file = Utils::outputPathBase(onSTBC) + dirs[i_dir] + "/particles_matched.root";
      std::unique_ptr<TFile> fParticles(TFile::Open(file));
      TTree* tp = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
      if (!tp) {
        std::cerr << "Warning: Cannot read " << file << ", skipping." << std::endl;
        continue;
      }

      std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr;
      std::vector<int> *pdg = nullptr, *nHits = nullptr;
      std::vector<unsigned>* generation = nullptr;
      std::vector<bool>* matched = nullptr;
      tp->SetBranchAddress("pt", &pt);
      tp->SetBranchAddress("eta", &eta);
      tp->SetBranchAddress("phi", &phi);
      tp->SetBranchAddress("particle_type", &pdg);
      tp->SetBranchAddress("number_of_hits", &nHits);
      tp->SetBranchAddress("generation", &generation);
      tp->SetBranchAddress("matched", &matched);
      for (Long64_t i = 0; i < tp->GetEntries(); i++) {
        tp->GetEntry(i);
        for (size_t k = 0; k < pt->size(); k++) {
          const auto it = std::find(pdgCodes.begin(), pdgCodes.end(), std::abs(pdg->at(k)));
          if (it == pdgCodes.end() || generation->at(k) != 0) continue;
          if (!acceptance.isAccepted(pt->at(k), eta->at(k), phi->at(k), std::abs(pdg->at(k)), nHits->at(k))) continue;
          const unsigned i_part = it - pdgCodes.begin();
          effPt[s][i_part]->Fill(matched->at(k), pt->at(k));
          if (pt->at(k) > ptCutForEta) effEta[s][i_part]->Fill(matched->at(k), std::abs(eta->at(k)));
        }
      }
    }

    // integrated efficiency above ptCutForEta, per species
    printf("%s, efficiency for p_T > %.2g GeV:", labels[s].Data(), ptCutForEta);
    for (unsigned i_part = 0; i_part < names.size(); i_part++) {
      const TH1* total = effEta[s][i_part]->GetTotalHistogram();
      const TH1* passed = effEta[s][i_part]->GetPassedHistogram();
      printf("  %s %.3f (%.0f)", names[i_part].Data(),
             total->GetEntries() ? passed->GetEntries() / total->GetEntries() : 0., total->GetEntries());
    }
    printf("\n");
  }

  // ---- plots: one pad per species, samples overlaid
  const TString base = Utils::aliceProjectsBase(onSTBC) + "acts-studies/geo_studies/";
  const TString figureDir = base + "figures/geo/" + tag + "/efficiency_forward/";
  gSystem->mkdir(figureDir, true);
  gSystem->mkdir(base + "histos/geo", true);
  auto draw = [&](std::vector<std::vector<TEfficiency*>>& effs, const TString& name, const TString& xTitle,
                  double xMin, double xMax, bool logx, const TString& selection) {
    TCanvas c("c_" + name, name, 1200, 900);
    c.Divide(2, 2);
    std::vector<TLegend*> legends;
    for (unsigned i_part = 0; i_part < names.size(); i_part++) {
      c.cd(i_part + 1);
      gPad->SetLogx(logx);
      gPad->SetGridy();
      gPad->DrawFrame(xMin, 0, xMax, 1.05)->SetTitle(symbols[i_part] + ", " + selection + ";" + xTitle + ";Efficiency");
      legends.push_back(new TLegend(0.55, 0.15, 0.88, 0.15 + 0.06 * samples.size()));
      legends.back()->SetBorderSize(0);
      legends.back()->SetFillStyle(0);
      for (unsigned s = 0; s < samples.size(); s++) {
        effs[s][i_part]->SetLineColor(colors[s % colors.size()]);
        effs[s][i_part]->SetMarkerColor(colors[s % colors.size()]);
        effs[s][i_part]->SetMarkerStyle(20);
        effs[s][i_part]->SetMarkerSize(0.8);
        effs[s][i_part]->Draw("P SAME");
        legends.back()->AddEntry(effs[s][i_part], labels[s], "lp");
      }
      legends.back()->Draw();
    }
    c.SaveAs(figureDir + name + ".pdf");
    for (TLegend* leg : legends) delete leg;
  };
  draw(effPt, "efficiency_vs_pT", "p_{T} (GeV/c)", 0.01, 10, true, Form("%.1f < |#eta| < %.1f", etaMin, etaMax));
  draw(effEta, "efficiency_vs_eta", "|#eta|", etaMin, etaMax, false, Form("p_{T} > %.2g GeV/c", ptCutForEta));

  TFile out(base + "histos/geo/eff_forward_" + tag + ".root", "RECREATE");
  for (unsigned s = 0; s < samples.size(); s++)
    for (unsigned i_part = 0; i_part < names.size(); i_part++) {
      effPt[s][i_part]->Write();
      effEta[s][i_part]->Write();
    }
  printf("Figures in %s, efficiencies in %s\n", figureDir.Data(), out.GetName());
}
