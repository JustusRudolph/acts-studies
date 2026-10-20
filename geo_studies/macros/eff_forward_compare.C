#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile.h>
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
 * pT > ptCutForEta. Denominator: particles in particles_matched.root (from
 * global/add_matching_to_particles.C) with at least minHits hits (7, as
 * alice3_reproduce.C). Numerator: those flagged matched.
 * Also the number of hits of these particles (all four species, pT > ptCutForEta): its
 * mean vs |eta| and its distribution.
 * Every figure is made separately for primaries (generation 0) and secondaries
 * (generation > 0).
 */

namespace EffForward {
const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<int> pdgCodes = {11, 211, 2212, 321};
const std::vector<TString> names = {"electron", "pion", "proton", "kaon"};
const std::vector<TString> symbols = {"e^{#pm}", "#pi^{#pm}", "p^{#pm}", "K^{#pm}"};
// indexed by isSecondary: generation 0, generation > 0
const std::vector<TString> origins = {"primaries", "secondaries"};
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

  // [isSecondary][sample][species]
  std::vector<std::vector<std::vector<TEfficiency*>>> effPt(origins.size(), std::vector<std::vector<TEfficiency*>>(samples.size()));
  auto effEta = effPt;
  // [isSecondary][sample]
  std::vector<std::vector<TProfile*>> hitsVsEta(origins.size());
  std::vector<std::vector<TH1D*>> hitsDist(origins.size());
  for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
    for (const bool isSecondary : {false, true}) {
      const char* org = origins[isSecondary].Data();  // index with the bool
      hitsVsEta[isSecondary].push_back(new TProfile(Form("hitsVsEta_%s_%u", org, i_sample), "", 20, etaMin, etaMax));
      hitsDist[isSecondary].push_back(new TH1D(Form("hitsDist_%s_%u", org, i_sample), "", 20, 0.5, 20.5));
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        effPt[isSecondary][i_sample].push_back(new TEfficiency(Form("effPt_%s_%s_%u", names[i_part].Data(), org, i_sample), "",
                                              nPtBins, ptBins.data()));
        effEta[isSecondary][i_sample].push_back(new TEfficiency(Form("effEta_%s_%s_%u", names[i_part].Data(), org, i_sample), "",
                                               20, etaMin, etaMax));
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
        printf("%s: %zu / %zu directories (%.0f%%)\n", labels[i_sample].Data(), i_dir, dirs.size(),
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
      for (Long64_t i_event = 0; i_event < tp->GetEntries(); i_event++) {
        tp->GetEntry(i_event);
        for (size_t i_particle = 0; i_particle < pt->size(); i_particle++) {
          const auto it = std::find(pdgCodes.begin(), pdgCodes.end(), std::abs(pdg->at(i_particle)));
          if (it == pdgCodes.end()) continue;
          if (!acceptance.isAccepted(pt->at(i_particle), eta->at(i_particle), phi->at(i_particle), std::abs(pdg->at(i_particle)), nHits->at(i_particle))) continue;
          const unsigned i_part = it - pdgCodes.begin();
          const bool isSecondary = generation->at(i_particle) > 0;
          effPt[isSecondary][i_sample][i_part]->Fill(matched->at(i_particle), pt->at(i_particle));
          if (pt->at(i_particle) > ptCutForEta) {
            effEta[isSecondary][i_sample][i_part]->Fill(matched->at(i_particle), std::abs(eta->at(i_particle)));
            hitsVsEta[isSecondary][i_sample]->Fill(std::abs(eta->at(i_particle)), nHits->at(i_particle));
            hitsDist[isSecondary][i_sample]->Fill(nHits->at(i_particle));
          }
        }
      }
    }

    // integrated efficiency above ptCutForEta, per species
    for (const bool isSecondary : {false, true}) {
      printf("%s, %s, efficiency for p_T > %.2g GeV:", labels[i_sample].Data(), origins[isSecondary].Data(), ptCutForEta);
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        const TH1* total = effEta[isSecondary][i_sample][i_part]->GetTotalHistogram();
        const TH1* passed = effEta[isSecondary][i_sample][i_part]->GetPassedHistogram();
        printf("  %s %.3f (%.0f)", names[i_part].Data(),
               total->GetEntries() ? passed->GetEntries() / total->GetEntries() : 0., total->GetEntries());
      }
      printf("\n");
    }
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
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        effs[i_sample][i_part]->SetLineColor(colors[i_sample % colors.size()]);
        effs[i_sample][i_part]->SetMarkerColor(colors[i_sample % colors.size()]);
        effs[i_sample][i_part]->SetMarkerStyle(20);
        effs[i_sample][i_part]->SetMarkerSize(0.8);
        effs[i_sample][i_part]->Draw("P SAME");
        legends.back()->AddEntry(effs[i_sample][i_part], labels[i_sample], "lp");
      }
      legends.back()->Draw();
    }
    c.SaveAs(figureDir + name + ".pdf");
    for (TLegend* leg : legends) delete leg;
  };
  auto style = [&](TH1* h, unsigned i_sample) {
    h->SetLineColor(colors[i_sample % colors.size()]);
    h->SetMarkerColor(colors[i_sample % colors.size()]);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
  };
  auto legend = [&](double x1, double y2, auto& hists) {
    TLegend* leg = new TLegend(x1, y2 - 0.06 * samples.size(), x1 + 0.25, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) leg->AddEntry(hists[i_sample], labels[i_sample], "lp");
    leg->Draw();
    return leg;
  };

  for (const bool isSecondary : {false, true}) {
    const TString& org = origins[isSecondary];
    draw(effPt[isSecondary], "efficiency_vs_pT_" + org, "p_{T} (GeV/c)", 0.01, 10, true,
         org + Form(", %.1f < |#eta| < %.1f", etaMin, etaMax));
    draw(effEta[isSecondary], "efficiency_vs_eta_" + org, "|#eta|", etaMin, etaMax, false,
         org + Form(", p_{T} > %.2g GeV/c", ptCutForEta));

    // ---- number of hits: mean vs |eta| with the difference to the first sample below, and the
    // distribution (unit area, including overflow), samples overlaid
    const TString hitsTitle = "e, #pi, p, K " + org +
                              Form(", %.1f < |#eta| < %.1f, p_{T} > %.2g GeV/c", etaMin, etaMax, ptCutForEta);
    {
      TCanvas c("c_nHits_vs_eta", "nHits vs eta", 800, 800);
      TPad top("top", "", 0, 0.3, 1, 1), bottom("bottom", "", 0, 0, 1, 0.3);
      top.SetBottomMargin(0.02);
      bottom.SetTopMargin(0.03);
      bottom.SetBottomMargin(0.3);
      top.SetGridy();
      bottom.SetGridy();
      top.Draw();
      bottom.Draw();

      top.cd();
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        hitsVsEta[isSecondary][i_sample]->SetTitle(hitsTitle + ";;Mean number of hits");
        style(hitsVsEta[isSecondary][i_sample], i_sample);
        hitsVsEta[isSecondary][i_sample]->SetMinimum(8);
        hitsVsEta[isSecondary][i_sample]->SetMaximum(15);
        hitsVsEta[isSecondary][i_sample]->GetXaxis()->SetLabelSize(0);
        hitsVsEta[isSecondary][i_sample]->Draw(i_sample == 0 ? "E" : "E SAME");
      }
      TLegend* leg = legend(0.65, 0.88, hitsVsEta[isSecondary]);

      // difference to the first sample, errors added as uncorrelated; bins empty in either
      // sample are left empty
      bottom.cd();
      TH1D* baseline = hitsVsEta[isSecondary][0]->ProjectionX("hitsVsEta_baseline_" + org);
      std::vector<TH1D*> diffs;
      double dMax = 0;
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        diffs.push_back(hitsVsEta[isSecondary][i_sample]->ProjectionX(Form("hitsVsEta_diff_%s_%u", org.Data(), i_sample)));
        for (int i_bin = 1; i_bin <= diffs[i_sample]->GetNbinsX(); i_bin++) {
          if (hitsVsEta[isSecondary][i_sample]->GetBinEntries(i_bin) == 0 || hitsVsEta[isSecondary][0]->GetBinEntries(i_bin) == 0) {
            diffs[i_sample]->SetBinContent(i_bin, 0);
            diffs[i_sample]->SetBinError(i_bin, 0);
            continue;
          }
          const double diff = diffs[i_sample]->GetBinContent(i_bin) - baseline->GetBinContent(i_bin);
          const double err = i_sample == 0 ? 0 : std::hypot(diffs[i_sample]->GetBinError(i_bin), baseline->GetBinError(i_bin));
          diffs[i_sample]->SetBinContent(i_bin, diff);
          diffs[i_sample]->SetBinError(i_bin, err);
          dMax = std::max(dMax, std::abs(diff) + err);
        }
      }
      const double dHalf = 1.1 * dMax + 1e-3;
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        style(diffs[i_sample], i_sample);
        diffs[i_sample]->SetTitle(";|#eta|;Diff. to " + labels[0]);
        diffs[i_sample]->SetMinimum(-dHalf);
        diffs[i_sample]->SetMaximum(dHalf);
        diffs[i_sample]->GetXaxis()->SetTitleSize(0.12);
        diffs[i_sample]->GetXaxis()->SetLabelSize(0.1);
        diffs[i_sample]->GetYaxis()->SetTitleSize(0.1);
        diffs[i_sample]->GetYaxis()->SetTitleOffset(0.45);
        diffs[i_sample]->GetYaxis()->SetLabelSize(0.09);
        diffs[i_sample]->GetYaxis()->SetNdivisions(505);
        diffs[i_sample]->Draw(i_sample == 0 ? "E" : "E SAME");
      }
      // the baseline itself is zero everywhere, drawn as a line in its colour
      TLine zero(diffs[0]->GetXaxis()->GetXmin(), 0, diffs[0]->GetXaxis()->GetXmax(), 0);
      zero.SetLineColor(colors[0]);
      zero.SetLineWidth(2);
      zero.Draw();
      c.SaveAs(figureDir + "nHits_vs_eta_" + org + ".pdf");
      delete leg;
    }
    {
      TCanvas c("c_nHits_distribution", "nHits distribution", 800, 600);
      gPad->SetGridy();
      double yMax = 0;
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        TH1D* h = hitsDist[isSecondary][i_sample];
        const double total = h->Integral(0, h->GetNbinsX() + 1);
        if (total > 0) h->Scale(1. / total);
        yMax = std::max(yMax, h->GetMaximum() + h->GetBinError(h->GetMaximumBin()));
      }
      for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
        hitsDist[isSecondary][i_sample]->SetTitle(hitsTitle + ";Number of hits;Fraction of particles");
        style(hitsDist[isSecondary][i_sample], i_sample);
        hitsDist[isSecondary][i_sample]->SetMinimum(0);
        hitsDist[isSecondary][i_sample]->SetMaximum(1.2 * yMax);
        hitsDist[isSecondary][i_sample]->Draw(i_sample == 0 ? "E" : "E SAME");
      }
      TLegend* leg = legend(0.15, 0.88, hitsDist[isSecondary]);
      c.SaveAs(figureDir + "nHits_distribution_" + org + ".pdf");
      delete leg;
    }
  }

  TFile out(base + "histos/geo/eff_forward_" + tag + ".root", "RECREATE");
  for (const bool isSecondary : {false, true})
    for (unsigned i_sample = 0; i_sample < samples.size(); i_sample++) {
      for (unsigned i_part = 0; i_part < names.size(); i_part++) {
        effPt[isSecondary][i_sample][i_part]->Write();
        effEta[isSecondary][i_sample][i_part]->Write();
      }
      hitsVsEta[isSecondary][i_sample]->Write();
      hitsDist[isSecondary][i_sample]->Write();
    }
  printf("Figures in %s, efficiencies in %s\n", figureDir.Data(), out.GetName());
}
