#include <iostream>
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>

#include "../../Utils/BasicFunctions.h"

// Distribution of the production vertex radius r = sqrt(vx^2 + vy^2) of all
// secondary (generation > 0) particles in <input_dir>/particles_simulation.root,
// input_dir being relative to the ACTSO2 output base
void secondary_vertex_r(const TString input_dir="pythia_2500ev/seed_0",
                        const bool onSTBC=false,
                        const double r_max=1500.,  // mm
                        const unsigned n_bins=300) {
  gStyle->SetOptStat(0);

  const TString sim_filename = Utils::outputPathBase(onSTBC) + input_dir + "/particles_simulation.root";
  TFile* sim_file = TFile::Open(sim_filename);
  if (!sim_file || sim_file->IsZombie()) {
    std::cout << "Error: Could not open " << sim_filename << std::endl;
    return;
  }
  TTree* sim_tree = (TTree*)sim_file->Get("particles");
  if (!sim_tree) {
    std::cout << "Error: particles tree not found in " << sim_filename << std::endl;
    sim_file->Close();
    return;
  }

  std::vector<float>* vx = nullptr;
  std::vector<float>* vy = nullptr;
  std::vector<unsigned>* generation = nullptr;
  sim_tree->SetBranchAddress("vx",         &vx);
  sim_tree->SetBranchAddress("vy",         &vy);
  sim_tree->SetBranchAddress("generation", &generation);

  TH1D* secondary_vertex_r_hist = new TH1D("secondary_vertex_r_hist",
                                           "Production vertex r of secondary particles (generation > 0);"
                                           "Vertex r (mm); Entries",
                                           n_bins, 0, r_max);

  unsigned long n_secondaries = 0;
  for (unsigned i_ev = 0; i_ev < sim_tree->GetEntries(); i_ev++) {
    sim_tree->GetEntry(i_ev);
    for (size_t i_mcp = 0; i_mcp < generation->size(); i_mcp++) {
      if (generation->at(i_mcp) == 0) continue;  // primaries
      secondary_vertex_r_hist->Fill(std::hypot(vx->at(i_mcp), vy->at(i_mcp)));
      n_secondaries++;
    }
  }  // loop over events
  std::cout << "Filled " << n_secondaries << " secondaries from "
            << sim_tree->GetEntries() << " events" << std::endl;

  // figures/<input_dir with "/" replaced by "_">
  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/pythia_studies/figures/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);

  TCanvas* canvas = new TCanvas("c_secondary_vertex_r", "Secondary vertex r", 800, 600);
  secondary_vertex_r_hist->SetLineColor(kBlue);
  secondary_vertex_r_hist->SetLineWidth(1);
  secondary_vertex_r_hist->GetXaxis()->SetTitleOffset(1.2);
  secondary_vertex_r_hist->Draw("HIST");
  gPad->SetLogy(1);
  canvas->SaveAs(Form("%s/secondary_vertex_r.pdf", output_dir.Data()));

  sim_file->Close();
}
