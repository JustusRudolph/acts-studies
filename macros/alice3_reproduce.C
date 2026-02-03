#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

#include "../utils/cuts.C"

void alice3_reproduce(const std::string input_data_type="default_seeds_2T_geo_oct25",
                      const std::string phi_eff_input_data_type="central_eta_001_1G_piplus_default_seeds",
                      const bool include_muons=false) {
  
  // Set ROOT style
  gStyle->SetOptStat(0);
  
  // Define particle types and pTs
  std::vector<TString> particles = {"electron", "piplus", "proton", "kaonplus"};
  std::vector<TString> particle_labels = {"e^{-}", "#pi^{+}", "p^{+}", "K^{+}"};
  std::vector<Int_t> colours_per_particle = {kMagenta, kBlack, kRed, kBlue};
  if (include_muons) {
    particles.push_back("muon");
    particle_labels.push_back("#mu^{-}");
    colours_per_particle.push_back(kGreen+3);
  }
  std::vector<int> pTs_mev = {50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                              2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
  // to have central bin values correct for plotting
  std::array<double, 26> pT_bin_edges_mev = {45., 55., 65., 75., 85., 95., 105., 195., 205., 395., 405.,
                                            595., 605., 795., 805., 995., 1005., 2995., 3005., 4995., 5005.,
                                            6995., 7005., 8995., 9005., 10995.};

  std::vector<TString> pT_strings = {"50M", "60M", "70M", "80M", "90M", "100M", "150M", "200M", "300M", "400M", "500M",
                                    "600M", "700M", "800M", "900M", "1G", "2G", "3G", "4G", "5G", "6G",
                                    "7G", "8G", "9G", "10G"};
  std::vector<TString> pT_labels = {"50 MeV/c", "60 MeV/c", "70 MeV/c", "80 MeV/c", "90 MeV/c", "100 MeV/c",
                                    "150 MeV/c", "200 MeV/c", "300 MeV/c", "400 MeV/c", "500 MeV/c", "600 MeV/c",
                                    "700 MeV/c", "800 MeV/c", "900 MeV/c", "1 GeV/c", "2 GeV/c", "3 GeV/c", "4 GeV/c",
                                    "5 GeV/c", "6 GeV/c", "7 GeV/c", "8 GeV/c", "9 GeV/c", "10 GeV/c"};

  std::vector<TEfficiency*> eff_hists_central(particles.size());  // |eta| < 1
  std::vector<TEfficiency*> eff_seeding_hists_central(particles.size());  // no cuts
  float abs_eta_cut = 1.0;                                          

  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    const TString particle = particles[i_part];
    eff_hists_central[i_part] = new TEfficiency(Form("eff_hist_%s", particle.Data()), Form("Efficiency vs p_{T} for %s", particle.Data()),
                                    pT_bin_edges_mev.size() - 1, pT_bin_edges_mev.data());
    eff_seeding_hists_central[i_part] = new TEfficiency(Form("eff_seeding_hist_%s", particle.Data()),
                                                        Form("Seeding Efficiency vs p_{T} for %s", particle.Data()),
                                                        pT_bin_edges_mev.size() - 1, pT_bin_edges_mev.data());

    for (unsigned i_pT = 0; i_pT < pT_strings.size(); i_pT++) {
      const TString energy = pT_strings[i_pT];
      TString dir_name = Form("data/%s/%s/%s", input_data_type.c_str(), particle.Data(), energy.Data());
      TString perf_ambi_filename = dir_name + "/performance_finding_ambi.root";
      TString perf_seed_filename = dir_name + "/performance_seeding.root";
      TFile* perf_ambi_file = TFile::Open(perf_ambi_filename);
      TFile* perf_seed_file = TFile::Open(perf_seed_filename);
      if (!perf_ambi_file || perf_ambi_file->IsZombie()) {
        std::cout << "Warning: Could not open " << perf_ambi_filename << std::endl;
        continue;
      }
      if (!perf_seed_file || perf_seed_file->IsZombie()) {
        std::cout << "Warning: Could not open " << perf_seed_filename << std::endl;
        continue;
      }

      TTree* ambi_mc_tree = (TTree*)perf_ambi_file->Get("matchingdetails");
      TTree* seed_mc_tree = (TTree*)perf_seed_file->Get("matchingdetails");

      std::vector<bool>* matched = nullptr;
      std::vector<bool>* isSecondary = nullptr;
      std::vector<float>* eta = nullptr;
      // std::vector<float>* phi = nullptr;

      std::vector<bool>* matched_seed = nullptr;
      std::vector<bool>* isSecondary_seed = nullptr;
      std::vector<float>* eta_seed = nullptr;

      ambi_mc_tree->SetBranchAddress("matched", &matched);
      ambi_mc_tree->SetBranchAddress("isSecondary", &isSecondary);
      ambi_mc_tree->SetBranchAddress("eta", &eta);
      // ambi_mc_tree->SetBranchAddress("phi", &phi);

      if (seed_mc_tree != nullptr) {
        seed_mc_tree->SetBranchAddress("matched", &matched_seed);
        seed_mc_tree->SetBranchAddress("isSecondary", &isSecondary_seed);
        seed_mc_tree->SetBranchAddress("eta", &eta_seed);
      }

      for (unsigned i_ev = 0; i_ev < ambi_mc_tree->GetEntries(); i_ev++) {
        ambi_mc_tree->GetEntry(i_ev);
        for (size_t i_mcp = 0; i_mcp < matched->size(); i_mcp++) {
          if (std::abs(eta->at(i_mcp)) < abs_eta_cut && !isSecondary->at(i_mcp)) {
            // TODO: UNCOMMENT THIS WHEN YOU HAVE PHI OUTPUT IN MATCHINGDETAILS
            // ALSO NOT RELEVANT FOR GEANT4 GEOMETRY (PERFECT CYLINDER)
            // if (Cuts::petalcut_phi(phi->at(i_mcp)))
            //   continue;  // if at petal, skip
            eff_hists_central[i_part]->Fill(matched->at(i_mcp), pTs_mev[i_pT]);
          }
        }
      }
      if (seed_mc_tree == nullptr)
        continue;  // no seeding efficiency available, so don't try to fill

      for (unsigned i_ev = 0; i_ev < seed_mc_tree->GetEntries(); i_ev++) {
        seed_mc_tree->GetEntry(i_ev);
        for (size_t i_mcp = 0; i_mcp < matched_seed->size(); i_mcp++) {
          if (std::abs(eta_seed->at(i_mcp)) < abs_eta_cut && !isSecondary_seed->at(i_mcp)) {
            eff_seeding_hists_central[i_part]->Fill(matched_seed->at(i_mcp), pTs_mev[i_pT]);
          }
        }
      }
    }  // loop over pT
  }  // loop over particles

  // ------------------- Now get efficiency vs phi for central eta -------------------
  TEfficiency* eff_wrt_phi_near_central = new TEfficiency("eff_wrt_phi_near_central", "Efficiency vs #phi;#phi;Efficiency",
                                             100, -3.14, 3.14);
  TEfficiency* eff_wrt_phi_complete_central = new TEfficiency("eff_wrt_phi_complete_central", "Efficiency vs #phi;#phi;Efficiency",
                                             100, -3.14, 3.14);
  auto near_central_eta_acceptance_phi_eff = [](float eta) {
    return (std::abs(eta) > 7e-3) && (std::abs(eta) < 1e-2);
  };
  auto complete_central_eta_acceptance_phi_eff = [](float eta) {
    return (std::abs(eta) < 1e-4);
  };
  TString phi_eff_path = Form("data/%s/performance_finding_ambi.root",
                              phi_eff_input_data_type.c_str());

  TFile* phi_eff_file = TFile::Open(phi_eff_path);
  if (!phi_eff_file || phi_eff_file->IsZombie()) {
    std::cout << "Warning: Could not open " << phi_eff_path << std::endl;
  }

  TTree* phi_eff_mctree = (TTree*)phi_eff_file->Get("matchingdetails");
  std::vector<bool>* matched_phi = nullptr;
  std::vector<float>* phi = nullptr;
  std::vector<float>* eta_phi = nullptr;
  phi_eff_mctree->SetBranchAddress("matched", &matched_phi);
  phi_eff_mctree->SetBranchAddress("phi", &phi);
  phi_eff_mctree->SetBranchAddress("eta", &eta_phi);

  for (unsigned i_ev = 0; i_ev < phi_eff_mctree->GetEntries(); i_ev++) {
    phi_eff_mctree->GetEntry(i_ev);
    for (size_t i_mcp = 0; i_mcp < matched_phi->size(); i_mcp++) {
      if (near_central_eta_acceptance_phi_eff(eta_phi->at(i_mcp))) {
        eff_wrt_phi_near_central->Fill(matched_phi->at(i_mcp), phi->at(i_mcp));
      }
      if (complete_central_eta_acceptance_phi_eff(eta_phi->at(i_mcp))) {
        eff_wrt_phi_complete_central->Fill(matched_phi->at(i_mcp), phi->at(i_mcp));
      }
    }
  }

  // Now plot efficiency vs pT for all particles
  TCanvas* canvas = new TCanvas("canvas_efficiency_alice3", "Efficiency vs pT", 1500, 1300);
  canvas->Divide(2,2);
  canvas->cd(1);
  
  // LEGEND
  // put legend in bottom right corner: same for ambi and seeding plots
  TLegend* leg_eff = new TLegend(0.6, 0.15, 0.85, 0.35 );
  leg_eff->SetBorderSize(0);
  leg_eff->SetFillStyle(0);
  
  // PLOTS - First draw to establish the frame
  eff_hists_central[0]->SetTitle("Efficiency vs p_{T} for various particles (|#eta| < 1); p_{T} (MeV/c); Efficiency");
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_hists_central[i_part]->SetMarkerColor(colours_per_particle[i_part]);
    eff_hists_central[i_part]->SetLineColor(colours_per_particle[i_part]);
    eff_hists_central[i_part]->SetMarkerStyle(20);
    eff_hists_central[i_part]->SetMarkerSize(1);
    eff_hists_central[i_part]->SetLineStyle(2);
    eff_hists_central[i_part]->SetLineWidth(1);
    leg_eff->AddEntry(eff_hists_central[i_part], particle_labels[i_part], "lp");
    if (i_part == 0) {
      eff_hists_central[i_part]->Draw("EXP");
    } else {
      eff_hists_central[i_part]->Draw("EXP SAME");
    }
  }
  gPad->SetLogx(1);
  gPad->Update();  // need to update before getting painted graph
  
  // BACKGROUND GRID LINES - Draw after establishing frame
  // Add vertical dotted black lines at specific momentum values
  std::vector<double> line_positions;
  
  // 50-100 MeV in steps of 10
  for (int pt = 50; pt <= 100; pt += 10) {
    line_positions.push_back(pt);
  }
  
  // 100-1000 MeV in steps of 100
  for (int pt = 200; pt <= 1000; pt += 100) {  // start from 200 to avoid duplicate at 100
    line_positions.push_back(pt);
  }
  
  // 1000-10000 MeV in steps of 1000
  for (int pt = 2000; pt <= 10000; pt += 1000) {  // start from 2000 to avoid duplicate at 1000
    line_positions.push_back(pt);
  }
  
  for (unsigned i_line = 0; i_line < line_positions.size(); i_line++) {
    double pt = line_positions[i_line];
    TLine* vline = new TLine(pt, 0, pt, 1.05);
    vline->SetLineStyle(kDotted);
    vline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    vline->SetLineWidth(1);
    vline->Draw("same");
  }

  // Add horizontal lines at steps of 0.2
  for (double eff = 0.2; eff <= 1.0; eff += 0.2) {
    TLine* hline = new TLine(45, eff, 11000, eff);
    hline->SetLineStyle(kDotted);
    hline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    hline->SetLineWidth(1);
    hline->Draw("same");
  }
  
  // Redraw efficiency plots on top of grid lines
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_hists_central[i_part]->Draw("EXP SAME");
  }
  
  leg_eff->Draw();

  eff_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);  // avoid axis title overwriting axis labels
  eff_hists_central[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  eff_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetRangeUser(45, 11000);

  // ----------- Now plot seeding efficiency vs pT for all particles -----------
  canvas->cd(2);
  eff_seeding_hists_central[0]->SetTitle("Seeding Efficiency vs p_{T} for various particles (|#eta| < 1); p_{T} (MeV/c); Seeding Efficiency");
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_seeding_hists_central[i_part]->SetMarkerColor(colours_per_particle[i_part]);
    eff_seeding_hists_central[i_part]->SetLineColor(colours_per_particle[i_part]);
    eff_seeding_hists_central[i_part]->SetMarkerStyle(20);
    eff_seeding_hists_central[i_part]->SetMarkerSize(1);
    eff_seeding_hists_central[i_part]->SetLineStyle(2);
    eff_seeding_hists_central[i_part]->SetLineWidth(1);
    if (i_part == 0) {
      eff_seeding_hists_central[i_part]->Draw("EXP");
    } else {
      eff_seeding_hists_central[i_part]->Draw("EXP SAME");
    }
  }
  gPad->SetLogx(1);
  gPad->Update();  // need to update before getting painted graph

  // straight lines in background now just like before
  for (unsigned i_line = 0; i_line < line_positions.size(); i_line++) {
    double pt = line_positions[i_line];
    TLine* vline = new TLine(pt, 0, pt, 1.05);
    vline->SetLineStyle(kDotted);
    vline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    vline->SetLineWidth(1);
    vline->Draw("same");
  }
  for (double eff = 0.2; eff <= 1.0; eff += 0.2) {
    TLine* hline = new TLine(45, eff, 11000, eff);
    hline->SetLineStyle(kDotted);
    hline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    hline->SetLineWidth(1);
    hline->Draw("same");
  }

  // Redraw efficiency plots on top of grid lines
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_seeding_hists_central[i_part]->Draw("EXP SAME");
  }

  leg_eff->Draw();

  eff_seeding_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);  // avoid axis title overwriting axis labels
  eff_seeding_hists_central[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  eff_seeding_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetRangeUser(45, 11000);

  // Now plot efficiency vs phi for central eta
  canvas->cd(3);
  eff_wrt_phi_near_central->SetMarkerColor(kViolet);
  eff_wrt_phi_near_central->SetLineColor(kViolet);
  eff_wrt_phi_near_central->SetMarkerStyle(20);
  eff_wrt_phi_near_central->SetMarkerSize(1);

  eff_wrt_phi_complete_central->SetMarkerColor(kOrange+7);
  eff_wrt_phi_complete_central->SetLineColor(kOrange+7);
  eff_wrt_phi_complete_central->SetMarkerStyle(24);
  eff_wrt_phi_complete_central->SetMarkerSize(1);
  eff_wrt_phi_near_central->SetTitle("Efficiency vs #phi for 1 GeV/c #pi^{#pm}; #phi (rad); Efficiency");
  eff_wrt_phi_near_central->Draw("EP");
  eff_wrt_phi_complete_central->Draw("EP SAME");

  gPad->Update();
  eff_wrt_phi_near_central->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);

  TLegend* leg2 = new TLegend(0.6, 0.25, 0.85, 0.35);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(eff_wrt_phi_near_central, "0.007 < |#eta| < 0.01", "lp");
  leg2->AddEntry(eff_wrt_phi_complete_central, "|#eta| < 10^{-4}", "lp");
  leg2->Draw();

  canvas->SaveAs(Form("figures/alice3_reproduced_plots/%s_%s.pdf",
                      input_data_type.c_str(), phi_eff_input_data_type.c_str()));
}