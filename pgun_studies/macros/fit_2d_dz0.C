#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

void fit_2d_dz0(const std::string particle_to_check="muon", unsigned energy_mev=1000,
                const std::string input_data_type="default_seeds_2T_geo_oct25") {
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);

  TString energy_str = energy_mev < 1000 ? Form("%dM", energy_mev) : Form("%dG", energy_mev / 1000);
  TString energy_label = energy_mev < 1000 ? Form("%d MeV", energy_mev) : Form("%d GeV", energy_mev / 1000);

  TString filename = Form("data/%s/%s/%s/performance_fitting_ambi.root",
                          input_data_type.c_str(), particle_to_check.c_str(), energy_str.Data());

  TFile* fit_file = TFile::Open(filename);
  if (!fit_file || fit_file->IsZombie()) {
    std::cout << "Warning: Could not open " << filename << std::endl;
    return;
  }
  // get histos
  TH2F* res_d0_vs_eta = (TH2F*) fit_file->Get("res_d0_vs_eta");
  TH2F* res_z0_vs_eta = (TH2F*) fit_file->Get("res_z0_vs_eta");

  // fit slices in eta at the boundary: 17th and 18th bin 1-indexed
  int bin_to_fit_outer = 17;
  int bin_to_fit_central = 18;

  // --------------------------- D0 fitting ---------------------------
  // fit the outer bin
  TF1* fit_func_outer_d0 = new TF1("fit_func_outer_d0", "gaus", -0.1, 0.1);
  TH1D* slice_outer_d0 = res_d0_vs_eta->ProjectionY("slice_outer_d0", bin_to_fit_outer, bin_to_fit_outer);
  slice_outer_d0->Scale(1.0 / slice_outer_d0->Integral());  // normalise
  slice_outer_d0->Fit("fit_func_outer_d0", "QN0");
  // fit the central bin
  TF1* fit_func_central_d0 = new TF1("fit_func_central_d0", "gaus", -0.1, 0.1);
  TH1D* slice_central_d0 = res_d0_vs_eta->ProjectionY("slice_central_d0", bin_to_fit_central, bin_to_fit_central);
  slice_central_d0->Scale(1.0 / slice_central_d0->Integral());  // normalise
  slice_central_d0->Fit("fit_func_central_d0", "QN0");

  // --------------------------- Z0 fitting ---------------------------
  // fit the outer bin
  TF1* fit_func_outer_z0 = new TF1("fit_func_outer_z0", "gaus", -0.1, 0.1);
  TH1D* slice_outer_z0 = res_z0_vs_eta->ProjectionY("slice_outer_z0", bin_to_fit_outer, bin_to_fit_outer);
  slice_outer_z0->Scale(1.0 / slice_outer_z0->Integral());  // normalise
  slice_outer_z0->Fit("fit_func_outer_z0", "QN0");
  // fit the central bin
  TF1* fit_func_central_z0 = new TF1("fit_func_central_z0", "gaus", -0.1, 0.1);
  TH1D* slice_central_z0 = res_z0_vs_eta->ProjectionY("slice_central_z0", bin_to_fit_central, bin_to_fit_central);
  slice_central_z0->Scale(1.0 / slice_central_z0->Integral());  // normalise
  slice_central_z0->Fit("fit_func_central_z0", "QN0");

  // plot results
  TCanvas* canvas = new TCanvas("canvas", "D0 and Z0 fits", 1200, 600);
  canvas->Divide(2,1);
  canvas->cd(1);
  slice_outer_d0->SetLineColorAlpha(kRed, 0.5);
  slice_outer_d0->SetTitle(Form("%s %s D0 residual fits; D0 residual [mm]; Entries (normalised)",
                        energy_label.Data(), particle_to_check.c_str()));
  slice_outer_d0->Draw("PE");
  slice_outer_d0->GetXaxis()->SetRangeUser(-0.1, 0.1);
  slice_outer_d0->GetYaxis()->SetRangeUser(0, std::max(slice_outer_d0->GetMaximum(), slice_central_d0->GetMaximum())*1.2);
  fit_func_outer_d0->SetLineColor(kRed);
  fit_func_outer_d0->SetLineStyle(2); // dashed
  fit_func_outer_d0->Draw("same");
  slice_central_d0->SetLineColorAlpha(kBlue, 0.5);
  slice_central_d0->Draw("PE same");
  fit_func_central_d0->SetLineColor(kBlue);
  fit_func_central_d0->SetLineStyle(2);
  fit_func_central_d0->Draw("same");
  TLegend* leg_d0 = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg_d0->SetBorderSize(0);
  leg_d0->SetFillStyle(0);
  leg_d0->AddEntry(slice_outer_d0, Form("Outer (#sigma = %.3f mm)", fit_func_outer_d0->GetParameter(2)), "lp");
  leg_d0->AddEntry(slice_central_d0, Form("Central (#sigma = %.3f mm)", fit_func_central_d0->GetParameter(2)), "lp");
  leg_d0->Draw();

  canvas->cd(2);
  slice_outer_z0->SetLineColorAlpha(kRed, 0.5);
  slice_outer_z0->SetTitle(Form("%s %s Z0 residual fits; Z0 residual [mm]; Entries (normalised)",
                        energy_label.Data(), particle_to_check.c_str()));
  slice_outer_z0->Draw("PE");
  slice_outer_z0->GetXaxis()->SetRangeUser(-0.1, 0.1);
  slice_outer_z0->GetYaxis()->SetRangeUser(0, std::max(slice_outer_z0->GetMaximum(), slice_central_z0->GetMaximum())*1.2);
  fit_func_outer_z0->SetLineColor(kRed);
  fit_func_outer_z0->SetLineStyle(2); // dashed
  fit_func_outer_z0->Draw("same");
  slice_central_z0->SetLineColorAlpha(kBlue, 0.5);
  slice_central_z0->Draw("PE same");
  fit_func_central_z0->SetLineColor(kBlue);
  fit_func_central_z0->SetLineStyle(2);
  fit_func_central_z0->Draw("same");
  TLegend* leg_z0 = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg_z0->SetBorderSize(0);
  leg_z0->SetFillStyle(0);
  leg_z0->AddEntry(slice_outer_z0, Form("Outer (#sigma = %.3f mm)", fit_func_outer_z0->GetParameter(2)), "lp");
  leg_z0->AddEntry(slice_central_z0, Form("Central (#sigma = %.3f mm)", fit_func_central_z0->GetParameter(2)), "lp");
  leg_z0->Draw();

  canvas->SaveAs(Form("figures/pgun_plots/fit_d0_z0_%s_%s_%s.pdf",
                      particle_to_check.c_str(), energy_str.Data(), input_data_type.c_str()));
}