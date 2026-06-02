#include <iostream>
#include <vector>
#include <map>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

#include "../utils/cuts.C"  // includes Constants.C

void phi_dep_eff(const std::vector<TString>& phi_eff_input_data_types,
                 TEfficiency* eff_wrt_phi_near_central,
                 TEfficiency* eff_wrt_phi_complete_central) {
  auto near_central_eta_acceptance_phi_eff = [](float eta) {
    return (std::abs(eta) > 7e-3) && (std::abs(eta) < 1e-2);
  };
  auto complete_central_eta_acceptance_phi_eff = [](float eta) {
    return (std::abs(eta) < 1e-4);
  };
  for (const auto& phi_eff_input_data_type : phi_eff_input_data_types) {
    TString phi_eff_path = Form("/%s/performance_finding_ambi.root",
                                phi_eff_input_data_type.c_str());

    TFile* phi_eff_file = TFile::Open(phi_eff_path);
    if (!phi_eff_file || phi_eff_file->IsZombie()) {
      std::cout << "Warning: Could not open " << phi_eff_path << std::endl;
      continue;
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
      }  // loop over mc particles in event
    }  // loop over entries (events) in file
    phi_eff_file->Close();
  }  // loop over files
}

void alice3_reproduce(const bool onSTBC=false,
                      const std::vector<TString> input_dirs={"50k_pythia_1"},
                      const std::vector<TString> phi_eff_input_dirs={},
                      const float abs_eta_max = 1.0,
                      const float pT_fraction_kept = 0.8,
                      const bool include_muons=false) {
  
  // Set ROOT style
  gStyle->SetOptStat(0);
  // need base of where ACTSO2 output is stored to access the input files
  const TString actso2_base = Constants::getACTSO2Base(onSTBC);
  // Define particle types and pTs
  std::vector<TString> particles = {"electron", "piplus", "proton", "kaonplus"};
  std::vector<TString> particle_labels = {"e^{-}", "#pi^{+}", "p^{+}", "K^{+}"};
  std::vector<int> pdg_codes = {11, 211, 2212, 321};
  std::map<int, unsigned> pdg_to_index;
  for (unsigned i = 0; i < pdg_codes.size(); i++) {
    pdg_to_index[pdg_codes[i]] = i;
  }
  std::vector<Int_t> colours_per_particle = {kMagenta, kBlack, kRed, kBlue};
  if (include_muons) {
    particles.push_back("muon");
    particle_labels.push_back("#mu^{-}");
    pdg_codes.push_back(13);
    pdg_to_index[13] = pdg_codes.size() - 1;
    colours_per_particle.push_back(kGreen+3);
  }
  // Create logarithmic bins from 0.01 GeV to 10 GeV (10 MeV to 10000 MeV)
  const int n_bins = 20;
  const double pt_min = 0.01;  // GeV
  const double pt_max = 10.0;
  const double log_min = std::log10(pt_min);
  const double log_max = std::log10(pt_max);
  const double log_step = (log_max - log_min) / n_bins;
  
  std::vector<double> pT_bin_edges_gev(n_bins + 1);
  for (int i = 0; i <= n_bins; i++) {
    pT_bin_edges_gev[i] = std::pow(10.0, log_min + i * log_step);
  }

  // set up efficiency histograms
  std::vector<TEfficiency*> eff_hists_central(particles.size());  // |eta| < 1
  std::vector<TEfficiency*> eff_seeding_hists_central(particles.size());  // no cuts

  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    const TString particle = particles[i_part];
    eff_hists_central[i_part] = new TEfficiency(Form("eff_hist_%s", particle.Data()), Form("Efficiency vs p_{T} for %s", particle.Data()),
                                    pT_bin_edges_gev.size() - 1, pT_bin_edges_gev.data());
    eff_seeding_hists_central[i_part] = new TEfficiency(Form("eff_seeding_hist_%s", particle.Data()),
                                                        Form("Seeding Efficiency vs p_{T} for %s", particle.Data()),
                                                        pT_bin_edges_gev.size() - 1, pT_bin_edges_gev.data());
  }

  for (const TString& input_dir : input_dirs) {
    TString dir_name = Form("%s/%s", actso2_base.Data(), input_dir.Data());
    TString sim_matched_filename = dir_name + "/particles_simulation_matched.root";
    TString perf_seed_filename = dir_name + "/performance_seeding.root";
    TFile* sim_matched_file = TFile::Open(sim_matched_filename);
    TFile* perf_seed_file = TFile::Open(perf_seed_filename);
    if (!sim_matched_file || sim_matched_file->IsZombie()) {
      std::cout << "Warning: Could not open " << sim_matched_filename << std::endl;
      continue;
    }
    if (!perf_seed_file || perf_seed_file->IsZombie()) {
      std::cout << "Warning: Could not open " << perf_seed_filename << std::endl;
      continue;
    }

    TTree* sim_tree = (TTree*)sim_matched_file->Get("particles");
    TTree* seed_mc_tree = (TTree*)perf_seed_file->Get("matchingdetails");

    std::vector<float>* eta = nullptr;
    std::vector<float>* pt = nullptr;
    std::vector<int>* pdg = nullptr;
    std::vector<unsigned>* generation = nullptr;
    std::vector<std::vector<unsigned>>* matchedIdxs = nullptr;

    std::vector<bool>* matched_seed = nullptr;
    std::vector<bool>* isSecondary_seed = nullptr;
    std::vector<float>* eta_seed = nullptr;
    std::vector<int>* pdg_seed = nullptr;
    std::vector<float>* pT_initial_seed = nullptr;
    std::vector<float>* pT_final_seed = nullptr;
    std::vector<float>* p_initial_seed = nullptr;
    std::vector<float>* p_final_seed = nullptr;

    sim_tree->SetBranchAddress("eta",                &eta);
    sim_tree->SetBranchAddress("pt",                 &pt);
    sim_tree->SetBranchAddress("<pdg_branch>",       &pdg);
    sim_tree->SetBranchAddress("generation",         &generation);
    sim_tree->SetBranchAddress("matched_track_idxs", &matchedIdxs);

    if (seed_mc_tree != nullptr) {
      seed_mc_tree->SetBranchAddress("matched", &matched_seed);
      seed_mc_tree->SetBranchAddress("isSecondary", &isSecondary_seed);
      seed_mc_tree->SetBranchAddress("eta", &eta_seed);
      seed_mc_tree->SetBranchAddress("pdg", &pdg_seed);
      seed_mc_tree->SetBranchAddress("pT_final", &pT_final_seed);
      seed_mc_tree->SetBranchAddress("p_final", &p_final_seed);
      seed_mc_tree->SetBranchAddress("pT_initial", &pT_initial_seed);
      seed_mc_tree->SetBranchAddress("p_initial", &p_initial_seed);
    }

    for (unsigned i_ev = 0; i_ev < sim_tree->GetEntries(); i_ev++) {
      sim_tree->GetEntry(i_ev);
      for (size_t i_mcp = 0; i_mcp < eta->size(); i_mcp++) {
        if (std::find(pdg_codes.begin(),
                      pdg_codes.end(),
                      pdg->at(i_mcp)) == pdg_codes.end())
          continue;

        unsigned i_part = pdg_to_index[pdg->at(i_mcp)];
        if (std::abs(eta->at(i_mcp)) < abs_eta_max
            && generation->at(i_mcp) == 0) {
          eff_hists_central[i_part]->Fill(!matchedIdxs->at(i_mcp).empty(), pt->at(i_mcp));
        }
      }
    }  // loop over sim_matched tree entries
    sim_matched_file->Close();

    if (seed_mc_tree != nullptr) {
      // seeding efficiency available
      for (unsigned i_ev = 0; i_ev < seed_mc_tree->GetEntries(); i_ev++) {
        seed_mc_tree->GetEntry(i_ev);
        for (size_t i_mcp = 0; i_mcp < matched_seed->size(); i_mcp++) {
          if (std::find(pdg_codes.begin(),
                        pdg_codes.end(),
                        pdg_seed->at(i_mcp)) == pdg_codes.end())
            continue;  // particle type not in our list, skip
          unsigned i_part = pdg_to_index[pdg_seed->at(i_mcp)];
          if (std::abs(eta_seed->at(i_mcp)) < abs_eta_max
              && !isSecondary_seed->at(i_mcp)
              && pT_final_seed->at(i_mcp) > pT_fraction_kept * pT_initial_seed->at(i_mcp)) {
            eff_seeding_hists_central[i_part]->Fill(matched_seed->at(i_mcp), pT_initial_seed->at(i_mcp));
          }
        }
      }  // loop over seeding tree entries
      perf_seed_file->Close();
    }  // if seeding tree exists
  }  // loop over input files

  // ------------------- Now get efficiency vs phi for central eta -------------------
  TEfficiency* eff_wrt_phi_near_central = new TEfficiency("eff_wrt_phi_near_central", "Efficiency vs #phi;#phi;Efficiency",
                                             100, -3.14, 3.14);
  TEfficiency* eff_wrt_phi_complete_central = new TEfficiency("eff_wrt_phi_complete_central", "Efficiency vs #phi;#phi;Efficiency",
                                             100, -3.14, 3.14);
  if (phi_eff_input_dirs.size() > 0) {
    phi_dep_eff(phi_eff_input_dirs, eff_wrt_phi_near_central, eff_wrt_phi_complete_central);
  } else {
    std::cout << "No phi efficiency input data type provided, skipping phi efficiency plot." << std::endl;
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
    eff_hists_central[i_part]->SetLineStyle(1);
    eff_hists_central[i_part]->SetLineWidth(1);
    leg_eff->AddEntry(eff_hists_central[i_part], particle_labels[i_part], "lp");
    if (i_part == 0) {
      eff_hists_central[i_part]->Draw("EP");
    } else {
      eff_hists_central[i_part]->Draw("EP SAME");
    }
  }
  gPad->SetLogx(1);
  gPad->Update();  // need to update before getting painted graph
  
  // BACKGROUND GRID LINES - Draw after establishing frame
  // Add vertical dotted black lines at powers of 10 and intermediate values
  std::vector<double> line_positions;
  
  // Add lines at 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10 GeV
  line_positions = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
  
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
    TLine* hline = new TLine(0.01, eff, 10.0, eff);
    hline->SetLineStyle(kDotted);
    hline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    hline->SetLineWidth(1);
    hline->Draw("same");
  }
  
  // Redraw efficiency plots on top of grid lines
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_hists_central[i_part]->Draw("EP SAME");
  }
  
  leg_eff->Draw();

  eff_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);  // avoid axis title overwriting axis labels
  eff_hists_central[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  eff_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetRangeUser(0.005, 11.0);

  // ----------- Now plot seeding efficiency vs pT for all particles -----------
  canvas->cd(2);
  eff_seeding_hists_central[0]->SetTitle("Seeding Efficiency vs p_{T} for various particles (|#eta| < 1); p_{T} (MeV/c); Seeding Efficiency");
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_seeding_hists_central[i_part]->SetMarkerColor(colours_per_particle[i_part]);
    eff_seeding_hists_central[i_part]->SetLineColor(colours_per_particle[i_part]);
    eff_seeding_hists_central[i_part]->SetMarkerStyle(20);
    eff_seeding_hists_central[i_part]->SetMarkerSize(1);
    eff_seeding_hists_central[i_part]->SetLineStyle(1);
    eff_seeding_hists_central[i_part]->SetLineWidth(1);
    if (i_part == 0) {
      eff_seeding_hists_central[i_part]->Draw("EP");
    } else {
      eff_seeding_hists_central[i_part]->Draw("EP SAME");
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
    TLine* hline = new TLine(0.01, eff, 10.0, eff);
    hline->SetLineStyle(kDotted);
    hline->SetLineColor(TColor::GetColorTransparent(kBlack, 0.4));  // 40% opacity
    hline->SetLineWidth(1);
    hline->Draw("same");
  }

  // Redraw efficiency plots on top of grid lines
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    eff_seeding_hists_central[i_part]->Draw("EP SAME");
  }

  leg_eff->Draw();

  eff_seeding_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);  // avoid axis title overwriting axis labels
  eff_seeding_hists_central[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  eff_seeding_hists_central[0]->GetPaintedGraph()->GetXaxis()->SetRangeUser(0.005, 11.0);

  // Now plot efficiency vs phi for central eta
  if (phi_eff_input_data_type != "") {
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
  }

  canvas->SaveAs(Form("figures/alice3_reproduced_plots/%s_%s.pdf",
                      input_data_type.c_str(), phi_eff_input_data_type.c_str()));
}