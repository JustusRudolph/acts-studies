#include <iostream>
#include <vector>
#include <sstream>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

void zoom_eta_binning(const std::string input_data_type="default_seeds_1T_geo_oct25",
                      const float abs_eta_max=1e-4,
                      const std::string particle="muon",
                      const unsigned pT_mev=100000) {
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);

  TString pT_str = pT_mev < 1000 ? Form("%dM", pT_mev) : Form("%dG", pT_mev / 1000);
  std::stringstream ss_scientific;
  ss_scientific << std::setprecision(1) << std::scientific << abs_eta_max;
  
  TString path_prefix = input_data_type;
  TString output_prefix = input_data_type + Form("%s_zoomed_eta", ss_scientific.str().c_str());
  if (particle.length()) {  // string is not empty
    path_prefix += Form("/%s/%s", particle.c_str(), pT_str.Data());
    output_prefix = Form("%s_%s", particle.c_str(), pT_str.Data()) + output_prefix;
  }
  
  TString filename = Form("data/%s/performance_finding_ambi.root", path_prefix.Data());
  TFile* perf_file = TFile::Open(filename);
  if (!perf_file || perf_file->IsZombie()) {
    std::cout << "Warning: Could not open " << filename << std::endl;
    return;
  }
  // everything is stored as a length 1 vector since we have only one particle per event
  std::vector<float>* eta = nullptr;
  std::vector<bool>* wasMatched = nullptr;
  
  TTree* perf_tree = (TTree*) perf_file->Get("matchingdetails");
  perf_tree->SetBranchAddress("eta", &eta);
  perf_tree->SetBranchAddress("matched", &wasMatched);
  

  // declare histogram(s)
  TEfficiency* eff_wrt_eta = new TEfficiency("eff_wrt_eta", "Efficiency vs #eta;#eta;Efficiency",
                                             50, -abs_eta_max, abs_eta_max);

  for (unsigned i_ev = 0; i_ev < perf_tree->GetEntries(); i_ev++) {
    perf_tree->GetEntry(i_ev);
    for (size_t i_mcp = 0; i_mcp < wasMatched->size(); i_mcp++) {
      if (std::abs(eta->at(i_mcp)) < abs_eta_max) {
        eff_wrt_eta->Fill(wasMatched->at(i_mcp), eta->at(i_mcp));
      }
    }
  }
  TCanvas* canvas = new TCanvas("canvas", "Efficiency vs #eta", 800, 600);
  eff_wrt_eta->Draw();

  canvas->SaveAs(Form("figures/eta_efficiency/%s.pdf", output_prefix.Data()));
  perf_file->Close();
}
