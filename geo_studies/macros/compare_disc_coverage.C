#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <vector>

/*
 * Overlays the hit rate vs r, the hit rate per stave and the hit efficiencies of
 * several histogram files made by disc_coverage_v2 (run that first for every sample).
 * Plots go to figures/geo/compare/{occupancy_layer,efficiency_layer}/.
 */

namespace CompareDiscCoverage {

constexpr int nDiscs = 6;
const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                 kOrange + 7, kCyan + 2, kViolet + 5};

TString label(int side, int disc) {
  return Form("%s_disc_%d", side == 0 ? "bwd" : "fwd", disc);
}
TString labelInPlot(int side, int disc) {
  return Form("%s disc %d (%s)", side == 0 ? "Backward" : "Forward", disc,
              disc < 3 ? "ML" : "OT");
}

template <class T>
std::vector<T*> getAll(const std::vector<TFile*>& files, const TString& name) {
  std::vector<T*> objs;
  for (unsigned f = 0; f < files.size(); f++) {
    T* obj = static_cast<T*>(files[f]->Get(name)->Clone(Form("%s_%u", name.Data(), f)));
    if (auto* h = dynamic_cast<TH1*>(obj)) h->SetDirectory(nullptr);
    objs.push_back(obj);
  }
  return objs;
}

TLegend* makeLegend(const std::vector<TObject*>& objs, const std::vector<TString>& labels) {
  TLegend* leg = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (unsigned f = 0; f < objs.size(); f++) leg->AddEntry(objs[f], labels[f], "lep");
  return leg;
}

// draw all histograms on the current pad, with a y range that fits all of them
void drawHists(std::vector<TH1D*> hists, const std::vector<TString>& labels,
               const TString& title, bool logy) {
  gPad->SetLeftMargin(0.15);  // y label cut off otherwise
  gPad->SetLogy(logy);
  double yMin = 0, yMax = 0;
  for (TH1D* h : hists) {
    yMax = std::max(yMax, h->GetMaximum());
    const double minPositive = h->GetMinimum(0);
    if (minPositive > 0 && (yMin == 0 || minPositive < yMin)) yMin = minPositive;
  }
  for (unsigned f = 0; f < hists.size(); f++) {
    hists[f]->SetTitle(title);
    hists[f]->SetLineColor(colors[f % colors.size()]);
    hists[f]->SetMarkerColor(colors[f % colors.size()]);
    hists[f]->SetMarkerStyle(20);
    hists[f]->SetMarkerSize(0.4);
    hists[f]->SetMinimum(logy ? 0.5 * yMin : 0);
    hists[f]->SetMaximum(logy ? 2 * yMax : 1.1 * yMax);
    hists[f]->Draw(f == 0 ? "EP" : "EP SAME");
  }
  makeLegend(std::vector<TObject*>(hists.begin(), hists.end()), labels)->Draw();
}

void drawEffs(std::vector<TEfficiency*> effs, const std::vector<TString>& labels,
              const TString& title) {
  gPad->SetLeftMargin(0.15);  // y label cut off otherwise
  for (unsigned f = 0; f < effs.size(); f++) {
    effs[f]->SetTitle(title);
    effs[f]->SetLineColor(colors[f % colors.size()]);
    effs[f]->SetMarkerColor(colors[f % colors.size()]);
    effs[f]->SetMarkerStyle(20);
    effs[f]->SetMarkerSize(0.6);
    effs[f]->Draw(f == 0 ? "AP" : "P SAME");
  }
  gPad->Update();  // the painted graph only exists after this
  effs[0]->GetPaintedGraph()->SetMinimum(0);
  effs[0]->GetPaintedGraph()->SetMaximum(1.05);
  makeLegend(std::vector<TObject*>(effs.begin(), effs.end()), labels)->Draw();
}

void save(TCanvas* c, const TString& path) {
  gSystem->mkdir(gSystem->DirName(path), true);
  c->SaveAs(path);
}

}  // namespace CompareDiscCoverage

/*
 * histFiles: histogram files of disc_coverage_v2, relative to histos/geo/ unless the
 * path is absolute. labels: legend entries, the file names without "disc_coverage_" if
 * empty. outTag is appended to every plot name so comparisons don't overwrite each other.
 */
void compare_disc_coverage(const std::vector<TString> histFiles,
                           std::vector<TString> labels = {},
                           const TString outTag = "",
                           const bool onSTBC = false) {
  using namespace CompareDiscCoverage;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  const TString base = onSTBC ? "/data/alice/jrudolph/" : "/home/justus/projects/";
  const TString geo_studies_base = base + "alice/acts-studies/geo_studies/";
  const TString occupancyDir = geo_studies_base + "figures/geo/compare/occupancy_layer/";
  const TString efficiencyDir = geo_studies_base + "figures/geo/compare/efficiency_layer/";
  const TString suffix = outTag.IsNull() ? TString(".pdf") : "_" + outTag + ".pdf";

  std::vector<TFile*> files;
  for (const TString& path : histFiles) {
    files.push_back(TFile::Open(path.BeginsWith("/") ? path
                                                     : geo_studies_base + "histos/geo/" + path));
  }
  if (labels.empty()) {
    for (const TString& path : histFiles) {
      TString name = gSystem->BaseName(path);
      name.ReplaceAll("disc_coverage_", "");
      name.ReplaceAll(".root", "");
      labels.push_back(name);
    }
  }

  // ---- hit rate per area vs r, one canvas per side
  for (int side = 0; side < 2; side++) {
    TCanvas* c = new TCanvas(Form("c_rate_%d", side), "Hit rate vs r", 1400, 900);
    c->Divide(3, 2);
    for (int disc = 0; disc < nDiscs; disc++) {
      c->cd(disc + 1);
      drawHists(getAll<TH1D>(files, "hHitRate_wrtR_" + label(side, disc)), labels,
                labelInPlot(side, disc) + ": hit rate vs r;r (mm);Hit rate (cm^{-2} s^{-1})",
                true);  // rate falls steeply with r
    }
    save(c, occupancyDir + "hitRate_wrt_R_" + (side == 0 ? "bwd" : "fwd") + suffix);
  }

  // ---- total hit rate per stave: one canvas each for ML and OT, bwd top, fwd bottom
  for (int det = 0; det < 2; det++) {  // 0 ML, 1 OT
    TCanvas* c = new TCanvas(Form("c_staveRate_%d", det), "Hit rate per stave", 1400, 900);
    c->Divide(3, 2);
    for (int side = 0; side < 2; side++) {
      for (int k = 0; k < 3; k++) {
        const int disc = 3 * det + k;
        c->cd(3 * side + k + 1);
        drawHists(getAll<TH1D>(files, "hStaveRate_" + label(side, disc)), labels,
                  labelInPlot(side, disc) + ": hit rate per stave;Stave ID;Hit rate (s^{-1})",
                  false);
      }
    }
    save(c, occupancyDir + "hitRate_per_stave_" + (det == 0 ? "ML" : "OT") + suffix);
  }

  // ---- hit efficiency by layer: nominal acceptance on top, tolerance ring below
  {
    TCanvas* c = new TCanvas("c_eff_layer", "Hit efficiency by layer", 1400, 1000);
    c->Divide(2, 2);
    const std::vector<TString> region = {"", "_outside_nominal_ring"};
    const std::vector<TString> regionInPlot = {"nominal acceptance", "tolerance region"};
    for (int r = 0; r < 2; r++) {
      for (int side = 0; side < 2; side++) {
        c->cd(2 * r + side + 1);
        const TString sideName = side == 0 ? "bwd" : "fwd";
        drawEffs(getAll<TEfficiency>(files, "effDisc_" + sideName + region[r]), labels,
                 TString(side == 0 ? "Backward" : "Forward") + " hit efficiency by layer, " +
                   regionInPlot[r] + ";Layer;Efficiency");
      }
    }
    save(c, efficiencyDir + "hit_efficiency_by_layer" + suffix);
  }

  // ---- hit efficiency vs r, one canvas per side
  for (int side = 0; side < 2; side++) {
    TCanvas* c = new TCanvas(Form("c_eff_r_%d", side), "Hit efficiency vs r", 1400, 900);
    c->Divide(3, 2);
    const TString sideName = side == 0 ? "bwd" : "fwd";
    for (int disc = 0; disc < nDiscs; disc++) {
      c->cd(disc + 1);
      drawEffs(getAll<TEfficiency>(files, "effDisc_" + sideName + "_" + label(side, disc)),
               labels, labelInPlot(side, disc) + ": efficiency vs r;r (mm);Efficiency");
    }
    save(c, efficiencyDir + "efficiency_vs_r_" + sideName + suffix);
  }

  for (TFile* f : files) f->Close();
}
