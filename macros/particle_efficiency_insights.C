#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>


void particle_efficiency_insights(const std::string input_data_type="default_seeds_1T_geo_oct25",
                                  const std::string particle="muon",
                                  const unsigned pT_mev=1000,
                                  const float max_eta_abs=4,
                                  const bool old_pid_scheme=false) {

  // Set ROOT style
  gStyle->SetOptStat(0);

  TString pT_str = pT_mev < 1000 ? Form("%dM", pT_mev) : Form("%dG", pT_mev / 1000);
  TString pT_label = pT_mev < 1000 ? Form("%d MeV/c", pT_mev) : Form("%d GeV/c", pT_mev / 1000);

  TString path_prefix = input_data_type;
  TString output_prefix = TString( Form("eta_%.1f_", max_eta_abs) ) + TString(input_data_type);
  if (particle.length()) {  // string is not empty
    path_prefix += Form("/%s/%s", particle.c_str(), pT_str.Data());
    output_prefix = Form("%s_%s_", particle.c_str(), pT_str.Data()) + output_prefix;
  }

  // load data
  TFile* perf_file = TFile::Open(Form("data/%s/performance_finding_ambi.root", path_prefix.Data()));
  // TFile* sim_file = TFile::Open(Form("data/%s/particles_simulation.root", path_prefix.Data()));

  if (!perf_file || perf_file->IsZombie()) {
    std::cout << "Warning: Could not open performance file for path_prefix " << path_prefix.Data() << std::endl;
    return;
  }

  TTree* perf_tree = (TTree*) perf_file->Get("matchingdetails");
  // TTree* sim_tree = (TTree*) sim_file->Get("particles");

  if (!perf_tree) {
    std::cout << "Warning: Could not find matchingdetails tree in performance file." << std::endl;
    return;
  }

  printf("Loaded performance file %s with %lld entries\n", 
         perf_file->GetName(), perf_tree->GetEntries());
  // Set up branches
  unsigned evNo;
  std::vector<unsigned long>* pid_perf = nullptr;
  std::vector<unsigned long>* pid_sim = nullptr;
  std::vector<bool>* wasMatched = nullptr;
  std::vector<std::vector<unsigned>>* matchedTrackIdxs = nullptr;
  std::vector<float>* eta = nullptr;
  std::vector<int>* nHits = nullptr;

  std::string pid_branch_name = old_pid_scheme ? "particle_id" : "particle_id_particle";
  perf_tree->SetBranchAddress("event_nr", &evNo);
  perf_tree->SetBranchAddress(pid_branch_name.c_str(), &pid_perf);
  perf_tree->SetBranchAddress("matched", &wasMatched);
  // perf_tree->SetBranchAddress("matchedTrackIdxs", &matchedTrackIdxs);
  perf_tree->SetBranchAddress("eta", &eta);
  perf_tree->SetBranchAddress("nHits", &nHits);

  // create eta range based on max and min eta values in the data
  float eta_min = -max_eta_abs;
  float eta_max = max_eta_abs;

  // create profiles: nHits wrt eta for matched & unmatched
  TProfile* prof_nHits_matched =
    new TProfile("prof_nHits_matched",
                 "nHits vs eta for matched particles; #eta; N_{Hits}",
                 50, eta_min, eta_max);
  TProfile* prof_nHits_unmatched =
    new TProfile("prof_nHits_unmatched",
                 "nHits vs eta for unmatched particles; #eta; N_{Hits}",
                 50, eta_min, eta_max);
  TProfile* prof_nHits_all =
    new TProfile("prof_nHits_all",
                 "nHits vs eta for all particles; #eta; N_{Hits}",
                 50, eta_min, eta_max);
  // Also create number of duplicates wrt eta
  TProfile* prof_nDuplicates =
    new TProfile("prof_nDuplicates",
                 "nDuplicates vs eta for matched particles; #eta; N_{Duplicates}",
                  50, eta_min, eta_max);

  TEfficiency* eff_wrt_nHits =
    new TEfficiency("nHitsMatched", "Efficiency wrt N_{Hits}; N_{Hits}; Entries",
                    20, -0.5, 19.5);

  // matched & unmatched eta distributions for different nHits
  TH1D* particle_7_hit_wrt_eta_matched =
    new TH1D("particle_7_hit_wrt_eta_matched",
              Form("%s 7-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_8_hit_wrt_eta_matched =
    new TH1D("particle_8_hit_wrt_eta_matched",
              Form("%s 8-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_9_hit_wrt_eta_matched =
    new TH1D("particle_9_hit_wrt_eta_matched",
              Form("%s 9-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_10_hit_wrt_eta_matched =
    new TH1D("particle_10_hit_wrt_eta_matched",
              Form("%s 10-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_11_hit_wrt_eta_matched =
    new TH1D("particle_11_hit_wrt_eta_matched",
              Form("%s 11-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_7_hit_wrt_eta_unmatched =
    new TH1D("particle_7_hit_wrt_eta_unmatched",
              Form("%s 7-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_8_hit_wrt_eta_unmatched =
    new TH1D("particle_8_hit_wrt_eta_unmatched",
              Form("%s 8-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_9_hit_wrt_eta_unmatched =
    new TH1D("particle_9_hit_wrt_eta_unmatched",
              Form("%s 9-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_10_hit_wrt_eta_unmatched =
    new TH1D("particle_10_hit_wrt_eta_unmatched",
              Form("%s 10-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);
  TH1D* particle_11_hit_wrt_eta_unmatched =
    new TH1D("particle_11_hit_wrt_eta_unmatched",
              Form("%s 11-hit #eta distribution; #eta; Frequency", particle.c_str()),
              50, eta_min, eta_max);

  TH1D* reco_hit_distribution =
    new TH1D("reco_hit_distribution",
              Form("%s reconstructed hit distribution; #eta; Frequency", particle.c_str()),
              20, -0.5, 19.5);
  TH1D* nonreco_hit_distribution =
    new TH1D("nonreco_hit_distribution",
              Form("%s non-reconstructed hit distribution; #eta; Frequency", particle.c_str()),
              20, -0.5, 19.5);

  // printf("Created profiles. Going through %lld entries\n", perf_tree->GetEntries());
  for (int i_ev = 0; i_ev < perf_tree->GetEntries(); i_ev++) {
    // printf("Processing entry %d. ", i_ev); 
    perf_tree->GetEntry(i_ev);
    // printf("Got data. ");
    unsigned nParticles = pid_perf->size();
    // printf("Have %u particles in this event. Go through.\n", nParticles);
    for (unsigned i_part = 0; i_part < nParticles; i_part++) {
      if (eta->at(i_part) < -max_eta_abs || eta->at(i_part) > max_eta_abs)
        continue;  // skip particles outside eta range

      if (wasMatched->at(i_part)) {
        prof_nHits_matched->Fill(eta->at(i_part), nHits->at(i_part));
        // prof_nDuplicates->Fill(eta->at(i_part), matchedTrackIdxs->at(i_part).size() - 1);
        reco_hit_distribution->Fill(nHits->at(i_part));
        // Now hits
        if (nHits->at(i_part) == 7) {
          particle_7_hit_wrt_eta_matched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 8) {
          particle_8_hit_wrt_eta_matched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 9) {
          particle_9_hit_wrt_eta_matched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 10) {
          particle_10_hit_wrt_eta_matched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 11) {
          particle_11_hit_wrt_eta_matched->Fill(eta->at(i_part));
        }
      } else {
        prof_nHits_unmatched->Fill(eta->at(i_part), nHits->at(i_part));
        nonreco_hit_distribution->Fill(nHits->at(i_part));
        // Now hits
        if (nHits->at(i_part) == 7) {
          particle_7_hit_wrt_eta_unmatched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 8) {
          particle_8_hit_wrt_eta_unmatched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 9) {
          particle_9_hit_wrt_eta_unmatched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 10) {
          particle_10_hit_wrt_eta_unmatched->Fill(eta->at(i_part));
        } else if (nHits->at(i_part) == 11) {
          particle_11_hit_wrt_eta_unmatched->Fill(eta->at(i_part));
        }
      }
      prof_nHits_all->Fill(eta->at(i_part), nHits->at(i_part));
      eff_wrt_nHits->Fill(wasMatched->at(i_part), nHits->at(i_part));

    }
  }
  // Plot everything
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1200, 1200);
  canvas->Divide(2,2);
  canvas->cd(1);
  gPad->SetLeftMargin(0.12);
  prof_nHits_matched->SetLineColor(kGreen+2);
  prof_nHits_unmatched->SetLineColor(kRed);
  prof_nHits_matched->SetMarkerColor(kGreen+2);
  prof_nHits_unmatched->SetMarkerColor(kRed);
  prof_nHits_all->SetLineColor(kBlack);
  prof_nHits_all->SetMarkerColor(kBlack);
  prof_nHits_matched->GetYaxis()->SetRangeUser(5., 16.);

  prof_nHits_matched->Draw("E1");
  prof_nHits_unmatched->Draw("E1 SAME");
  prof_nHits_all->Draw("E1 SAME");

  TLegend* leg1 = new TLegend(0.4, 0.7, 0.7, 0.82 );
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(prof_nHits_matched, "Reconstructed", "lp");
  leg1->AddEntry(prof_nHits_unmatched, "Non-reconstructed", "lp");
  leg1->AddEntry(prof_nHits_all, "All", "lp");
  leg1->Draw();

  canvas->cd(2);
  gPad->SetLeftMargin(0.12);
  eff_wrt_nHits->SetLineColor(kBlue);
  eff_wrt_nHits->SetMarkerColor(kBlue);
  eff_wrt_nHits->Draw("AP");

  canvas->cd(3);
  gPad->SetLeftMargin(0.12);
  TH1D* full_hit_distribution =
    (TH1D*) reco_hit_distribution->Clone("full_hit_distribution");
  full_hit_distribution->Add(nonreco_hit_distribution);
  full_hit_distribution->SetLineColor(kBlack);
  full_hit_distribution->SetMarkerColor(kBlack);
  reco_hit_distribution->SetLineColor(kGreen+2);
  reco_hit_distribution->SetMarkerColor(kGreen+2);
  nonreco_hit_distribution->SetLineColor(kRed);
  nonreco_hit_distribution->SetMarkerColor(kRed);
  full_hit_distribution->Draw("HIST");
  reco_hit_distribution->Draw("HIST SAME");
  nonreco_hit_distribution->Draw("HIST SAME");
  full_hit_distribution->SetTitle(Form("%s hit distribution; N_{Hits}; Frequency", particle.c_str()));
  full_hit_distribution->GetXaxis()->SetRangeUser(5, 15);
  gPad->SetLogy(1);

  TLegend* leg_hit_distr = new TLegend(0.2, 0.7, 0.4, 0.85);
  leg_hit_distr->SetBorderSize(0);
  leg_hit_distr->SetFillStyle(0);
  leg_hit_distr->AddEntry(full_hit_distribution, "All reconstructed", "l");
  leg_hit_distr->AddEntry(reco_hit_distribution, "Reconstructed", "l");
  leg_hit_distr->AddEntry(nonreco_hit_distribution, "Non-reconstructed", "l");
  leg_hit_distr->Draw();

  canvas->cd(4);
  gPad->SetLeftMargin(0.15);  // Increase left margin for y-axis label
  // For now want to draw just 8 & 11 hits, but change here as you want
  particle_8_hit_wrt_eta_matched->Scale(1. / particle_8_hit_wrt_eta_matched->Integral());
  particle_8_hit_wrt_eta_unmatched->Scale(1. / particle_8_hit_wrt_eta_unmatched->Integral());
  particle_11_hit_wrt_eta_matched->Scale(1. / particle_11_hit_wrt_eta_matched->Integral());
  particle_11_hit_wrt_eta_unmatched->Scale(1. / particle_11_hit_wrt_eta_unmatched->Integral());

  particle_8_hit_wrt_eta_matched->SetMarkerColor(kGreen+2);
  particle_8_hit_wrt_eta_matched->SetLineColor(kGreen+2);
  particle_8_hit_wrt_eta_unmatched->SetMarkerColor(kRed);
  particle_8_hit_wrt_eta_unmatched->SetLineColor(kRed);
  particle_11_hit_wrt_eta_matched->SetMarkerColor(kGreen+2);
  particle_11_hit_wrt_eta_matched->SetLineColor(kGreen+2);
  particle_11_hit_wrt_eta_unmatched->SetMarkerColor(kRed);
  particle_11_hit_wrt_eta_unmatched->SetLineColor(kRed);

  particle_8_hit_wrt_eta_matched->SetMarkerStyle(24);
  particle_8_hit_wrt_eta_unmatched->SetMarkerStyle(24);
  particle_11_hit_wrt_eta_matched->SetMarkerStyle(22);
  particle_11_hit_wrt_eta_unmatched->SetMarkerStyle(22);

  particle_11_hit_wrt_eta_matched->SetTitle(Form("%s #eta distribution for particles with various N_{Hits}; #eta; Density",
                                                particle.c_str()));
  particle_11_hit_wrt_eta_matched->Draw("PE");
  particle_11_hit_wrt_eta_unmatched->Draw("PE SAME");
  // particle_8_hit_wrt_eta_matched->Draw("PE SAME");
  particle_8_hit_wrt_eta_unmatched->Draw("PE SAME");
  particle_11_hit_wrt_eta_matched->GetYaxis()->SetRangeUser(0., 0.07);

  TLegend* leg_nHit_matched_eta_distr = new TLegend(0.16, 0.7, 0.5, 0.8);
  leg_nHit_matched_eta_distr->SetBorderSize(0);
  leg_nHit_matched_eta_distr->SetFillStyle(0);
  // leg_nHit_matched_eta_distr->AddEntry(particle_8_hit_wrt_eta_matched, "8 hits - Matched", "p");
  leg_nHit_matched_eta_distr->AddEntry(particle_8_hit_wrt_eta_unmatched, "8 hits - Unmatched", "p");
  leg_nHit_matched_eta_distr->AddEntry(particle_11_hit_wrt_eta_matched, "11 hits - Matched", "p");
  leg_nHit_matched_eta_distr->AddEntry(particle_11_hit_wrt_eta_unmatched, "11 hits - Unmatched", "p");
  leg_nHit_matched_eta_distr->Draw();

  canvas->SaveAs(Form("figures/particle_efficiency_insights/%s.pdf", output_prefix.Data()));
  delete canvas;

  // write histograms too
  printf("Writing output root file\n");
  TFile* out_file = TFile::Open(Form("hists/particle_efficiency_insights/%s.root", output_prefix.Data()), "RECREATE");
  prof_nHits_matched->Write();
  prof_nHits_unmatched->Write();
  prof_nDuplicates->Write();
  out_file->Close();
  delete out_file;

  perf_file->Close();
  // sim_file->Close();
  delete perf_file;
  // delete sim_file;
}