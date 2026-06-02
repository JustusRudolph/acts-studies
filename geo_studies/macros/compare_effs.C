#include <iostream>
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>

// ---- configure histograms to compare here ----
const std::vector<TString> kFindingHistNames = {
  "trackeff_vs_eta",
  "trackeff_vs_phi",
};

const std::vector<TString> kFindingTitles = {
  "Comparison in Track efficiency vs #eta",
  "Comparison in Track efficiency vs #phi",
};

const std::vector<TString> kFindingXLabels = {
  "#eta",
  "#phi (rad)",
};

const std::vector<TString> kFindingYLabels = {
  "Track efficiency",
  "Track efficiency",
};

const std::vector<TString> kFittingHistNames = {
  "reswidth_d0_vs_eta",
  "res_qopt",
  "reswidth_qopt_vs_eta",
};

const std::vector<TString> kFittingTitles = {
  "Comparison in d_{0} resolution vs #eta",
  "Comparison in q/p_{T} resolution vs q/p_{T}",
  "Comparison in q/p_{T} resolution vs #eta",
};

const std::vector<TString> kFittingHistXLabels = {
  "#eta",
  "q/p_{T} (1/GeV)",
  "#eta",
};

const std::vector<TString> kFittingHistYLabels = {
  "d_{0} Resolution (1/GeV)",
  "Frequency",
  "q/p_{T} Resolution (1/GeV)",
};
// ----------------------------------------------

TGraphErrors* makeRatioGraph(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {
  int n = std::min(g1->GetN(), g2->GetN());
  TGraphErrors* ratio = new TGraphErrors(n);
  for (int i = 0; i < n; i++) {
    double x, y1, y2, dummy;
    g1->GetPoint(i, x, y1);
    g2->GetPoint(i, dummy, y2);
    double r  = (y1 > 0.) ? y2 / y1 : 0.;
    double e1 = (g1->GetErrorYhigh(i) + g1->GetErrorYlow(i)) / 2.;
    double e2 = (g2->GetErrorYhigh(i) + g2->GetErrorYlow(i)) / 2.;
    double er = (y1 > 0. && y2 > 0.) ? r * std::sqrt((e1/y1)*(e1/y1) + (e2/y2)*(e2/y2)) : 0.;
    ratio->SetPoint(i, x, r);
    ratio->SetPointError(i, 0., er);
  }
  return ratio;
}

TH1* makeRatioHist(TH1* h1, TH1* h2) {
  TH1* ratio = (TH1*)h2->Clone("ratio");
  ratio->Divide(h1);
  return ratio;
}

void styleTopPad(TVirtualPad* pad) {
  pad->SetBottomMargin(0.03);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
}

void styleBottomPad(TVirtualPad* pad) {
  pad->SetTopMargin(0.03);
  pad->SetBottomMargin(0.35);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
}

void styleRatioAxis(TAxis* xaxis, TAxis* yaxis, const TString& xtitle) {
  double scale = 1. / 0.3;
  xaxis->SetTitle(xtitle);
  xaxis->SetTitleSize(0.04 * scale);
  xaxis->SetLabelSize(0.035 * scale);
  xaxis->SetTitleOffset(0.9);
  yaxis->SetTitle("Ratio");
  yaxis->SetTitleSize(0.04 * scale);
  yaxis->SetLabelSize(0.035 * scale);
  yaxis->SetTitleOffset(0.45);
  yaxis->SetNdivisions(504);
  yaxis->SetRangeUser(0.5, 1.5);
}

void drawHistograms(
  TFile* f1, TFile* f2,
  const std::vector<TString>& histNames,
  const std::vector<TString>& titles,
  const std::vector<TString>& xLabels,
  const std::vector<TString>& yLabels,
  TCanvas* c,
  const TString& outname,
  const TString& label1,
  const TString& label2
) {
  for (size_t idx = 0; idx < histNames.size(); ++idx) {
    const TString& hname  = histNames[idx];
    const TString& title  = titles[idx];
    const TString& xlabel = xLabels[idx];
    const TString& ylabel = yLabels[idx];
    TObject* obj1 = f1->Get(hname);
    TObject* obj2 = f2->Get(hname);

    if (!obj1) { std::cerr << "Not found in file 1: " << hname << std::endl; continue; }
    if (!obj2) { std::cerr << "Not found in file 2: " << hname << std::endl; continue; }

    c->Clear();
    TPad* padTop = new TPad("padTop", "", 0., 0.3, 1., 1.);
    TPad* padBot = new TPad("padBot", "", 0., 0.,  1., 0.3);
    styleTopPad(padTop);
    styleBottomPad(padBot);
    padTop->Draw();
    padBot->Draw();

    // ---- TH2 branch ----
    if (obj1->InheritsFrom(TH2::Class())) {
      padTop->Delete(); padBot->Delete();
      c->Clear();
      c->SetRightMargin(0.18);
      c->SetLeftMargin(0.12);
      c->SetBottomMargin(0.12);

      TH2* h1_2d = (TH2*)((TH2*)obj1)->Clone("h1_2d");
      TH2* h2_2d = (TH2*)((TH2*)obj2)->Clone("h2_2d");
      TH2* ratio2d = (TH2*)h2_2d->Clone("ratio_2d");
      ratio2d->Divide(h1_2d);
      ratio2d->SetTitle(Form("%s  ratio (%s / %s)", title.Data(), label2.Data(), label1.Data()));
      ratio2d->GetXaxis()->SetTitle(xlabel);
      ratio2d->GetYaxis()->SetTitle(ylabel);
      ratio2d->Draw("COLZ");

    // ---- TEfficiency branch ----
    } else if (obj1->InheritsFrom(TEfficiency::Class())) {
      TEfficiency* eff1 = (TEfficiency*)((TEfficiency*)obj1)->Clone("eff1");
      TEfficiency* eff2 = (TEfficiency*)((TEfficiency*)obj2)->Clone("eff2");

      eff1->SetLineColor(kBlue);  eff1->SetMarkerColor(kBlue);  eff1->SetMarkerStyle(20);
      eff2->SetLineColor(kRed);   eff2->SetMarkerColor(kRed);   eff2->SetMarkerStyle(21);

      padTop->cd();
      eff1->SetTitle(Form("%s; ; %s", title.Data(), ylabel.Data()));
      eff1->Draw("AP");
      gPad->Update();
      eff1->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);
      eff1->GetPaintedGraph()->GetXaxis()->SetLabelSize(0.);
      eff2->Draw("P SAME");
      gPad->Update();

      TLegend* leg = new TLegend(0.55, 0.15, 0.88, 0.30);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry(eff1, label1.Data(), "lp");
      leg->AddEntry(eff2, label2.Data(), "lp");
      leg->Draw();
      gPad->Update();

      TGraphAsymmErrors* g1 = (TGraphAsymmErrors*)eff1->GetPaintedGraph()->Clone();
      TGraphAsymmErrors* g2 = (TGraphAsymmErrors*)eff2->GetPaintedGraph()->Clone();

      padBot->cd();
      TGraphErrors* ratio = makeRatioGraph(g1, g2);
      ratio->SetMarkerStyle(20); ratio->SetMarkerColor(kBlack); ratio->SetLineColor(kBlack);
      ratio->SetTitle(Form("; %s; Ratio", xlabel.Data()));
      ratio->Draw("AP");
      gPad->Update();
      styleRatioAxis(ratio->GetXaxis(), ratio->GetYaxis(), xlabel);
      gPad->Update();

      double xmin = g1->GetXaxis()->GetXmin();
      double xmax = g1->GetXaxis()->GetXmax();
      TLine* line = new TLine(xmin, 1., xmax, 1.);
      line->SetLineStyle(2); line->SetLineColor(kGray+1);
      line->Draw();

    // ---- TH1 branch ----
    } else if (obj1->InheritsFrom(TH1::Class())) {
      TH1* h1 = (TH1*)((TH1*)obj1)->Clone("h1");
      TH1* h2 = (TH1*)((TH1*)obj2)->Clone("h2");

      h1->SetLineColor(kBlue);  h1->SetMarkerColor(kBlue);  h1->SetMarkerStyle(20);
      h2->SetLineColor(kRed);   h2->SetMarkerColor(kRed);   h2->SetMarkerStyle(21);

      padTop->cd();
      h1->SetTitle(Form("%s; ; %s", title.Data(), ylabel.Data()));
      h1->GetXaxis()->SetLabelSize(0.);
      double ymax = std::max(h1->GetMaximum(), h2->GetMaximum()) * 1.15;
      h1->GetYaxis()->SetRangeUser(0., ymax);
      h1->Draw("EP");
      h2->Draw("EP SAME");

      TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry(h1, label1.Data(), "lp");
      leg->AddEntry(h2, label2.Data(), "lp");
      leg->Draw();
      gPad->Update();

      padBot->cd();
      TH1* ratio = makeRatioHist(h1, h2);
      ratio->SetLineColor(kBlack); ratio->SetMarkerColor(kBlack); ratio->SetMarkerStyle(20);
      ratio->SetTitle(Form("; %s; Ratio", xlabel.Data()));
      styleRatioAxis(ratio->GetXaxis(), ratio->GetYaxis(), xlabel);
      ratio->Draw("EP");
      gPad->Update();

      double xmin = ratio->GetXaxis()->GetXmin();
      double xmax = ratio->GetXaxis()->GetXmax();
      TLine* line = new TLine(xmin, 1., xmax, 1.);
      line->SetLineStyle(2); line->SetLineColor(kGray+1);
      line->Draw();

    } else {
      std::cerr << "Unsupported type for " << hname << ": " << obj1->ClassName() << std::endl;
      continue;
    }

    c->Print(outname);
  }
}

void compare_effs(
  const TString path1,
  const TString path2,
  const TString label1 = "Sample 1",
  const TString label2 = "Sample 2"
) {
  gStyle->SetOptStat(0);

  const TString base = "/home/justus/projects/alice/ACTSO2/output/";

  TFile* fFind1 = TFile::Open(base + path1 + "/performance_finding_ambi.root");
  TFile* fFind2 = TFile::Open(base + path2 + "/performance_finding_ambi.root");
  TFile* fFit1  = TFile::Open(base + path1 + "/performance_fitting_ambi.root");
  TFile* fFit2  = TFile::Open(base + path2 + "/performance_fitting_ambi.root");

  if (!fFind1 || fFind1->IsZombie()) { std::cerr << "Cannot open finding file 1" << std::endl; return; }
  if (!fFind2 || fFind2->IsZombie()) { std::cerr << "Cannot open finding file 2" << std::endl; return; }
  if (!fFit1  || fFit1->IsZombie())  { std::cerr << "Cannot open fitting file 1" << std::endl; return; }
  if (!fFit2  || fFit2->IsZombie())  { std::cerr << "Cannot open fitting file 2" << std::endl; return; }

  TString outname = Form("figures/compare/compare_%s_vs_%s.pdf", label1.Data(), label2.Data());

  TCanvas* c = new TCanvas("c", "comparison", 800, 700);
  c->Print(outname + "[");

  drawHistograms(fFind1, fFind2, kFindingHistNames, kFindingTitles, kFindingXLabels, kFindingYLabels, c, outname, label1, label2);
  drawHistograms(fFit1,  fFit2,  kFittingHistNames, kFittingTitles, kFittingHistXLabels, kFittingHistYLabels, c, outname, label1, label2);

  c->Print(outname + "]");
  std::cout << "Saved: " << outname << std::endl;

  fFind1->Close(); fFind2->Close();
  fFit1->Close();  fFit2->Close();
}
