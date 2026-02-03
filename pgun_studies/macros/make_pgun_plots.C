#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

void make_pgun_plots(const std::string particle_to_check="muon",
                     const std::string input_data_type="truth_seeds_2T_geo_oct25") {
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);
  std::vector<int> energies_mev = {100, 1000, 10000, 100000};
  std::vector<TString> energies = {"100M", "1G", "10G", "100G"};
  std::vector<TString> energyLabels = {"100 MeV", "1 GeV", "10 GeV", "100 GeV"};
  std::vector<Int_t> colours_per_energy = {kRed, kBlue, kGreen+2, kYellow+1};
  std::vector<TProfile*> n_hits_wrt_eta(energies.size());
  std::vector<TProfile*> e_loss_wrt_eta(energies.size());
  
  for (size_t i = 0; i < energies.size(); i++) {
    n_hits_wrt_eta[i] = new TProfile(Form("n_hits_wrt_eta_%s", energies[i].Data()),
                                     Form("Number of hits vs #eta for %s;#eta;N_{hits}", energyLabels[i].Data()),
                                     30, -4, 4);
    e_loss_wrt_eta[i] = new TProfile(Form("e_loss_wrt_eta_%s", energies[i].Data()),
                                     Form("Energy loss vs #eta for %s;#eta;Energy loss [MeV]", energyLabels[i].Data()),
                                     30, -4, 4);
    TString filename = Form("data/%s/%s/%s/particles_simulation.root",
                            input_data_type.c_str(), particle_to_check.c_str(), energies[i].Data());
    TFile* mcfile = TFile::Open(filename);
    if (!mcfile || mcfile->IsZombie()) {
      std::cout << "Warning: Could not open " << filename << std::endl;
      continue;
    }
    std::cout << "Opened file: " << filename << std::endl;
    // everything is stored as a length 1 vector since we have only one particle per event
    std::vector<float>* eta = nullptr;
    std::vector<unsigned>* nHits = nullptr;
    std::vector<float>* eLoss = nullptr;
    TTree* mc_tree = (TTree*)mcfile->Get("particles");
    mc_tree->SetBranchAddress("number_of_hits", &nHits);
    mc_tree->SetBranchAddress("e_loss", &eLoss);
    mc_tree->SetBranchAddress("eta", &eta);
    std::cout << "Set branch addresses" << std::endl;

    unsigned nEntries = mc_tree->GetEntries();
    unsigned nEtas{0}, nNHits{0}, nELoss{0};
    for (unsigned i_mcEv = 0; i_mcEv < nEntries; i_mcEv++) {
      // printf("Checking entry %d\n", i_mcEv);
      mc_tree->GetEntry(i_mcEv);
      // printf("Filling entry %d\n", i_mcEv);
      for (size_t i_mcp = 0; i_mcp < eta->size(); i_mcp++) {
        n_hits_wrt_eta[i]->Fill(eta->at(i_mcp), nHits->at(i_mcp));
        e_loss_wrt_eta[i]->Fill(eta->at(i_mcp), eLoss->at(i_mcp));
      }
      nEtas+=eta->size();
      nNHits+=nHits->size();
      nELoss+=eLoss->size();
    }
    std::cout << "Number of etas: " << nEtas << std::endl;
    std::cout << "Number of nHits: " << nNHits << std::endl;
    std::cout << "Number of eLoss: " << nELoss << std::endl;
    mcfile->Close();
    delete mcfile;
  }
  std::cout << "Filled n_hits_wrt_eta histograms" << std::endl;
  // Now time for d0 and z0 resolution plots
  std::vector<TH1F*> d0_resolution(energies.size());
  std::vector<TH1F*> z0_resolution(energies.size());
  for (size_t i_energy = 0; i_energy < energies.size(); i_energy++) {
    TString filename = Form("data/%s/%s/%s/performance_fitting_ambi.root",
                            input_data_type.c_str(), particle_to_check.c_str(), energies[i_energy].Data());
    TFile* perf_fitting_file = TFile::Open(filename);
    if (!perf_fitting_file || perf_fitting_file->IsZombie()) {
      std::cout << "Warning: Could not open " << filename << std::endl;
      continue;
    }
    TH1F* d0_resolution_hist = (TH1F*)perf_fitting_file->Get("reswidth_d0_vs_eta");
    TH1F* z0_resolution_hist = (TH1F*)perf_fitting_file->Get("reswidth_z0_vs_eta");
    d0_resolution[i_energy] = (TH1F*)d0_resolution_hist->Clone(Form("d0_resolution_%s", energies[i_energy].Data()));
    d0_resolution[i_energy]->SetTitle(Form("D0 resolution vs #eta for %s", energyLabels[i_energy].Data()));
    d0_resolution[i_energy]->SetYTitle("D0 resolution [mm]");
    d0_resolution[i_energy]->SetDirectory(0);
    z0_resolution[i_energy] = (TH1F*)z0_resolution_hist->Clone(Form("z0_resolution_%s", energies[i_energy].Data()));
    z0_resolution[i_energy]->SetTitle(Form("Z0 resolution vs #eta for %s", energyLabels[i_energy].Data()));
    z0_resolution[i_energy]->SetYTitle("Z0 resolution [mm]");
    z0_resolution[i_energy]->SetDirectory(0);
    printf("Cloned d0 and z0 resolution histograms.\n");
    perf_fitting_file->Close();
    delete perf_fitting_file;
  }
  std::cout << "Filled d0 and z0 resolution histograms" << std::endl;
  std::cout << "Writing histograms to file" << std::endl;
  TFile* outputFile = new TFile(Form("hists/pgun_plots_%s.root", particle_to_check.c_str()), "recreate");
  std::cout << "Opened output file" << std::endl;
  for (unsigned i_energy = 0; i_energy < energies.size(); i_energy++) {
    printf("Writing histograms for energy %s\n", energies[i_energy].Data());
    n_hits_wrt_eta[i_energy]->Write();
    z0_resolution[i_energy]->Write();
    d0_resolution[i_energy]->Write();
  }
  outputFile->Close();
  delete outputFile;

  // Now make the plots in addition to the hists that are written separately
  TLegend* leg_res = new TLegend(0.425, 0.6, 0.625, 0.75);
  leg_res->SetBorderSize(0);
  leg_res->SetFillStyle(0);

  for (unsigned i_energy = 0; i_energy < energies.size(); i_energy++) {
    // set colours of markers and lines
    n_hits_wrt_eta[i_energy]->SetMarkerColor(colours_per_energy[i_energy]);
    n_hits_wrt_eta[i_energy]->SetLineColor(colours_per_energy[i_energy]);
    e_loss_wrt_eta[i_energy]->SetMarkerColor(colours_per_energy[i_energy]);
    e_loss_wrt_eta[i_energy]->SetLineColor(colours_per_energy[i_energy]);
    z0_resolution[i_energy]->SetMarkerColor(colours_per_energy[i_energy]);
    z0_resolution[i_energy]->SetLineColor(colours_per_energy[i_energy]);
    d0_resolution[i_energy]->SetMarkerColor(colours_per_energy[i_energy]);
    d0_resolution[i_energy]->SetLineColor(colours_per_energy[i_energy]);
    leg_res->AddEntry(d0_resolution[i_energy], energyLabels[i_energy].Data(), "p");
  }
  TCanvas* canvas = new TCanvas("canvas", "PGUN Plots", 1200, 800);
  canvas->Divide(2, 2);
  canvas->cd(1);
  n_hits_wrt_eta[0]->Draw("P");
  n_hits_wrt_eta[0]->GetYaxis()->SetRangeUser(7, 15);
  for (unsigned i_energy = 0; i_energy < energies.size(); i_energy++) {
    n_hits_wrt_eta[i_energy]->Draw("P SAME");
  }
  leg_res->Draw();
  // energy loss
  canvas->cd(2);
  e_loss_wrt_eta[0]->SetTitle(Form("Energy loss vs #eta for %s;#eta;Energy loss [MeV]",
                                   particle_to_check.c_str()));
  e_loss_wrt_eta[0]->Draw("P");
  e_loss_wrt_eta[0]->GetYaxis()->SetRangeUser(0., 0.5);
  for (unsigned i_energy = 0; i_energy < energies.size(); i_energy++) {
    e_loss_wrt_eta[i_energy]->Draw("P SAME");
  }
  leg_res->Draw();
  
  // resolution plots
  canvas->cd(3);
  d0_resolution[0]->SetTitle(Form("D_{0} resolution vs #eta for %s;#eta;D_{0} resolution [mm]",
                                  particle_to_check.c_str()));
  d0_resolution[0]->Draw("P");
  d0_resolution[0]->GetYaxis()->SetRangeUser(0.001, 0.5);
  for (unsigned i_energy = 1; i_energy < energies.size(); i_energy++) {
    d0_resolution[i_energy]->Draw("P SAME");
  }
  TLegend* leg_res_higher_up = (TLegend*) leg_res->Clone("leg_res_higher_up");
  leg_res_higher_up->SetY1(0.7);
  leg_res_higher_up->SetY2(0.85);
  leg_res_higher_up->Draw();
  gPad->SetLogy(1);
  canvas->cd(4);
  z0_resolution[0]->Draw("P");
  z0_resolution[0]->GetYaxis()->SetRangeUser(0.001, 2);
  z0_resolution[0]->SetTitle(Form("Z_{0} resolution vs #eta for %s;#eta;Z_{0} resolution [mm]",
                                  particle_to_check.c_str()));
  for (unsigned i_energy = 1; i_energy < energies.size(); i_energy++) {
    z0_resolution[i_energy]->Draw("P SAME");
  }
  leg_res->Draw();
  gPad->SetLogy(1);
  // save the canvas
  canvas->SaveAs(Form("figures/pgun_plots/%s_%s.pdf", particle_to_check.c_str(), input_data_type.c_str()));
}