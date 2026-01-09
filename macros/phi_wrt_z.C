/*
 * This macro is a very simple one, only intended for single particle runs. It
 * computes the azimuthal angle phi of all hit positions and plots it wrt z.
 */
#include <iostream>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

void phi_wrt_z(const TString particle="piplus", const unsigned pT_mev=10000) {

  // Set ROOT style, no auto-legend
  gStyle->SetOptStat(0);

  TString pT_str = pT_mev < 1000 ? Form("%dM", pT_mev) : Form("%dG", pT_mev / 1000);
  TString pT_label = pT_mev < 1000 ? Form("%d MeV/c", pT_mev) : Form("%d GeV/c", pT_mev / 1000);

  // load data, assume we are in pgun_studies
  TFile* hit_file = TFile::Open("../output_tutorial/reco_output_gun/hits.root");

  TTree* hit_tree = (TTree*) hit_file->Get("hits");

  printf("Loaded hit file %s with %lld entries\n", 
         hit_file->GetName(), hit_tree->GetEntries());

  // Set up branches
  float x, y, z;
  int hit_index;
  hit_tree->SetBranchAddress("tx", &x);
  hit_tree->SetBranchAddress("ty", &y);
  hit_tree->SetBranchAddress("tz", &z);
  hit_tree->SetBranchAddress("index", &hit_index);

  TGraph* phi_wrt_idx_graph = new TGraph();
  phi_wrt_idx_graph->SetTitle(Form("%s hit #phi vs hit index; hit index; #phi [rad]", particle.Data()));

  TGraph* z_wrt_idx_graph = new TGraph();

  unsigned nEntries = hit_tree->GetEntries();
  for (unsigned i_hit = 0; i_hit < nEntries; i_hit++) {
    hit_tree->GetEntry(i_hit);
    float phi = atan2f(y, x);
    phi_wrt_idx_graph->SetPoint(i_hit, hit_index, phi);
    z_wrt_idx_graph->SetPoint(i_hit, hit_index, z);
  }

  // Plot
  TCanvas* canvas = new TCanvas("canvas", "Phi and Z vs Hit Index", 1200, 800);
  canvas->SetLeftMargin(0.15);   // Make room for left y-axis label
  canvas->SetRightMargin(0.15);  // Make room for right axis
  
  // Draw phi graph first
  phi_wrt_idx_graph->SetMarkerStyle(20);
  phi_wrt_idx_graph->SetMarkerColor(kBlue);
  phi_wrt_idx_graph->Draw("AP");
  
  // Get the data ranges
  Double_t xmin, xmax, ymin_phi, ymax_phi;
  phi_wrt_idx_graph->ComputeRange(xmin, ymin_phi, xmax, ymax_phi);
  
  Double_t ymin_z, ymax_z;  
  z_wrt_idx_graph->ComputeRange(xmin, ymin_z, xmax, ymax_z);
  
  // Scale z data to phi coordinate system
  TGraph* z_scaled = new TGraph();
  for (int i = 0; i < z_wrt_idx_graph->GetN(); i++) {
    Double_t x_val, z_val;
    z_wrt_idx_graph->GetPoint(i, x_val, z_val);
    // Linear scaling: map z range to phi range
    Double_t z_mapped = ymin_phi + (z_val - ymin_z) * (ymax_phi - ymin_phi) / (ymax_z - ymin_z);
    z_scaled->SetPoint(i, x_val, z_mapped);
  }
  
  // Draw scaled z data
  z_scaled->SetMarkerStyle(21);
  z_scaled->SetMarkerColor(kRed);
  z_scaled->Draw("P same");
  
  // get pad margins
  Double_t leftMargin = canvas->GetLeftMargin();
  Double_t rightMargin = canvas->GetRightMargin();
  Double_t bottomMargin = canvas->GetBottomMargin();
  Double_t topMargin = canvas->GetTopMargin();
  // Create overlay pad for right axis  
  TPad* overlay = new TPad("overlay", "", 0, 0, 1, 1);
  overlay->SetFillStyle(4000);  // Transparent
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();
  
  // Calculate frame coordinates in pad normalised coordinates
  Double_t rightX = 1.0 - rightMargin;   // Right edge of plot frame
  Double_t bottomY = bottomMargin;       // Bottom edge of plot frame
  Double_t topY = 1.0 - topMargin;       // Top edge of plot frame
  
  // Draw TGaxis at the exact right frame edge
  TGaxis* rightaxis = new TGaxis(rightX, bottomY, rightX, topY,
                                ymin_z, ymax_z, 510, "+L");
  rightaxis->SetTitle("z [mm]");
  rightaxis->SetTitleColor(kRed);
  rightaxis->SetLabelColor(kRed);
  rightaxis->SetLineColor(kRed);
  rightaxis->SetTitleOffset(1.2);
  rightaxis->Draw();

  
  canvas->SaveAs(Form("figures/central_eta_studies/phi_wrt_idx_%s_%s.pdf",
                 particle.Data(), pT_str.Data()));
}