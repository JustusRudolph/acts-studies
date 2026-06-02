#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

#include "../utils/constants.C"

void plot_bremsstrahlung(const TString input_dir = "test_nEv20000_PID11_nMeasMin7_ckfChi2Meas45_ckfMeasPerSurf1",
                         const float max_eta_abs=1.0) {

  gStyle->SetOptStat(0);  // disable stats box
  // This file assumes that you will only look in the files that are given in the directory
  TString full_path_to_dir = Paths::data_base_dir + input_dir;

  TFile* particles_file = TFile::Open(full_path_to_dir + "/particles_simulation.root");
  if (!particles_file || particles_file->IsZombie()) {
    std::cout << "Warning: Could not open " << full_path_to_dir + "/particles_simulation.root" << std::endl;
    return;
  }

  // Now get ambi performance file to get matching details
  TFile* perf_ambi_file = TFile::Open(full_path_to_dir + "/performance_finding_ambi.root");
  if (!perf_ambi_file || perf_ambi_file->IsZombie()) {
    std::cout << "Warning: Could not open " << full_path_to_dir + "/performance_finding_ambi.root" << std::endl;
    return;
  }

  TTree* particles_tree = (TTree*)particles_file->Get("particles");
  TTree* perf_ambi_tree = (TTree*)perf_ambi_file->Get("matchingdetails");

  std::vector<int>* pdg_sim = nullptr;
  std::vector<int>* pdg_perf = nullptr;
  std::vector<float>* eta_sim = nullptr;
  std::vector<float>* eta_perf = nullptr;
  std::vector<float>* pT_initial = nullptr;
  std::vector<float>* pT_final = nullptr;
  std::vector<float>* pT_sim = nullptr;  // in particles_simulation.root
  std::vector<float>* p_initial = nullptr;
  std::vector<float>* p_final = nullptr;
  std::vector<float>* p_sim = nullptr;  // in particles_simulation.root
  std::vector<unsigned>* generation_sim = nullptr;
  std::vector<unsigned>* generation_perf = nullptr;
  std::vector<float>* v_x_sim = nullptr;
  std::vector<float>* v_y_sim = nullptr;
  std::vector<float>* v_z_sim = nullptr;

  // create log pt bins for gamma spectra
  const unsigned nBins = 100;
  double pTlogBins[nBins + 1];
  double ratio_pTlogBins[nBins + 1];
  double pt_min = 0.1;  // 10 MeV/c
  double pt_max = 10.0;  // 10 GeV/c
  double pTlogMin = std::log10(pt_min);
  double pTlogMax = std::log10(pt_max);
  for (unsigned i = 0; i <= nBins; ++i) {
    pTlogBins[i] = std::pow(10, pTlogMin + i * (pTlogMax - pTlogMin) / nBins);
    ratio_pTlogBins[i] = pTlogBins[i] / pt_max;
  }

  // create histogram of electron if they experienced bremsstrahlung or not, and if yes, how much
  // pT they lost in the first generation
  // Also create a histogram of total number of bremsstrahlung photons emitted per primary electron

  TH1D* total_pT_loss_hist = new TH1D("total_pT_loss_hist",
                                     Form("Total ratio p_{T} loss for primary electrons at |#eta| < %.1f;\
                                           p_{T} loss / e^{-} p_{T}; Fraction of electrons", max_eta_abs),
                                     100, ratio_pTlogBins);
  TH1D* brem_pT_loss_hist = new TH1D("brem_pT_loss_hist",
                                     Form("Relative p_{T} loss due to all bremsstrahlung for primary electrons at |#eta| < %.1f;\
                                           p_{T} loss / e^{-} p_{T}; Fraction of electrons", max_eta_abs),
                                     100, ratio_pTlogBins);
  TH1D* brem_pT_loss_hist_gen0 = new TH1D("brem_pT_loss_hist_gen0",
                                     Form("Relative p_{T} loss due to first bremsstrahlung for primary electrons at |#eta| < %.1f;\
                                           p_{T} loss / e^{-} p_{T}; Fraction of electrons", max_eta_abs),
                                     100, ratio_pTlogBins);
  TH1D* brem_photon_count_hist = new TH1D("brem_photon_count_hist",
                                         Form("Number of bremsstrahlung photons emitted per primary electron at |#eta| < %.1f;\
                                               Number of photons; Fraction of electrons", max_eta_abs),
                                         10, -0.5, 9.5);
  TEfficiency* brem_prob_wrt_pT_hist = new TEfficiency("brem_prob_wrt_pT_hist",
                                                  Form("Probability of bremsstrahlung emission for primary electrons at |#eta| < %.1f;\
                                                        p_{T} (GeV/c); Probability", max_eta_abs),
                                                        100, pTlogBins);
  TProfile* nPhotons_wrt_pT_profile = new TProfile("nPhotons_wrt_pT_profile",
                                                  Form("Average number of bremsstrahlung photons emitted for primary electrons at |#eta| < %.1f;\
                                                        p_{T} (GeV/c); Average number of photons", max_eta_abs),
                                                        100, pTlogBins);
  TH2D* gen1_photon_vertex_rz_hist = new TH2D("gen1_photon_vertex_rz_hist",
                                     Form("Position of first generation bremsstrahlung photon vertices for primary electrons at |#eta| < %.1f;\
                                           Vertex z (mm); Vertex r (mm)", max_eta_abs),
                                     200, -2500, 2500, 200, 0, 1000);
  TH2D* gen1_photon_vertex_rz_hist_zoom = new TH2D("gen1_photon_vertex_rz_hist_zoom",
                                     Form("Position of first generation bremsstrahlung photon vertices (zoomed) for primary electrons at |#eta| < %.1f;\
                                           Vertex z (mm); Vertex r (mm)", max_eta_abs),
                                     100, -100, 100, 100, 0, 100);
  TH2D* gen1_photon_vertex_rz_hist_highPtLoss = new TH2D("gen1_photon_vertex_rz_hist_highPtLoss",
                                     Form("Position of first generation bremsstrahlung photon vertices (p_{T} loss #geq 50%%) for primary electrons at |#eta| < %.1f;\
                                           Vertex z (mm); Vertex r (mm)", max_eta_abs),
                                     200, -2500, 2500, 200, 0, 1000);
  TH2D* gen1_photon_vertex_rz_hist_medPtLoss = new TH2D("gen1_photon_vertex_rz_hist_medPtLoss",
                                     Form("Position of first generation bremsstrahlung photon vertices (p_{T} loss 10-50%%) for primary electrons at |#eta| < %.1f;\
                                           Vertex z (mm); Vertex r (mm)", max_eta_abs),
                                     200, -2500, 2500, 200, 0, 1000);  
  particles_tree->SetBranchAddress("particle_type", &pdg_sim);
  particles_tree->SetBranchAddress("eta", &eta_sim);
  particles_tree->SetBranchAddress("vx", &v_x_sim);
  particles_tree->SetBranchAddress("vy", &v_y_sim);
  particles_tree->SetBranchAddress("vz", &v_z_sim);
  particles_tree->SetBranchAddress("pt", &pT_sim);
  particles_tree->SetBranchAddress("p", &p_sim);
  particles_tree->SetBranchAddress("generation", &generation_sim);

  perf_ambi_tree->SetBranchAddress("pdg", &pdg_perf);
  perf_ambi_tree->SetBranchAddress("eta", &eta_perf);
  perf_ambi_tree->SetBranchAddress("pT_initial", &pT_initial);
  perf_ambi_tree->SetBranchAddress("pT_final", &pT_final);
  perf_ambi_tree->SetBranchAddress("p_initial", &p_initial);
  perf_ambi_tree->SetBranchAddress("p_final", &p_final);
  perf_ambi_tree->SetBranchAddress("particle_id_generation", &generation_perf);

  for (unsigned i_ev = 0; i_ev < particles_tree->GetEntries(); i_ev++) {
    printf("Processing event %u / %lld: ", i_ev + 1, particles_tree->GetEntries());
    particles_tree->GetEntry(i_ev);
    perf_ambi_tree->GetEntry(i_ev);

    bool first_brem_found = false;
    unsigned n_brems = 0;  // to count number of photons from primary electron
    float total_brem_pT_loss = 0;  // to sum up total pT loss from bremsstrahlung for primary electron
    float primary_electron_pT_initial = 0;
    TString to_print = "";
    printf("Event has %lu simulated particles and %lu reconstructible particle(s)\n",
           pdg_sim->size(), pdg_perf->size());
    for (unsigned i_part = 0; i_part < pdg_perf->size(); i_part++) {
      if (std::abs(eta_perf->at(i_part)) > max_eta_abs)
      continue;  // skip particles outside eta range
      
      if (pdg_perf->at(i_part) == 11 && generation_perf->at(i_part) == 0) {  // primary electron
        float pT_loss = pT_initial->at(i_part) - pT_final->at(i_part);
        primary_electron_pT_initial = pT_initial->at(i_part);
        if (primary_electron_pT_initial > 0) {
          float relative_pT_loss = pT_loss / primary_electron_pT_initial;
          total_pT_loss_hist->Fill(relative_pT_loss);
        }
      }
    }
    
    // Go through simulated particles to fill brem histograms
    for (unsigned i_part = 0; i_part < pdg_sim->size(); i_part++) {
      to_print += Form("\tParticle %u out of %lu...", i_part + 1, pdg_sim->size());
      if (std::abs(eta_sim->at(i_part)) > max_eta_abs)
        continue;  // skip particles outside eta range

      to_print += " and is a " + TString(pdg_sim->at(i_part) == 11 ? "electron" : (pdg_sim->at(i_part) == 22 ? "photon" : "other particle"));
      if (pdg_sim->at(i_part) == 11 && generation_sim->at(i_part) == 0) {  // primary electron
        primary_electron_pT_initial = pT_sim->at(i_part);
      }
      // first generation bremsstrahlung photon with at least 1 MeV pT
      if (pdg_sim->at(i_part) == 22 && generation_sim->at(i_part) == 1 && pT_sim->at(i_part) > 0.001) {
        float relative_pT_loss = 0.0;
        if (primary_electron_pT_initial > 0) {
          relative_pT_loss = pT_sim->at(i_part) / primary_electron_pT_initial;
        }
        
        if (!first_brem_found) {
          // pT of the photon is the the pT loss of the electron
          if (primary_electron_pT_initial > 0) {
            brem_pT_loss_hist_gen0->Fill(relative_pT_loss);
          }
          first_brem_found = true;
        }
        total_brem_pT_loss += pT_sim->at(i_part);
        float vertex_r = std::sqrt(v_x_sim->at(i_part) * v_x_sim->at(i_part) +
                                   v_y_sim->at(i_part) * v_y_sim->at(i_part));
        float vertex_z = v_z_sim->at(i_part);
        gen1_photon_vertex_rz_hist->Fill(vertex_z, vertex_r);
        gen1_photon_vertex_rz_hist_zoom->Fill(vertex_z, vertex_r);
        
        // Fill pT loss category histograms
        if (relative_pT_loss >= 0.5) {
          gen1_photon_vertex_rz_hist_highPtLoss->Fill(vertex_z, vertex_r);
        } else if (relative_pT_loss >= 0.1 && relative_pT_loss < 0.5) {
          gen1_photon_vertex_rz_hist_medPtLoss->Fill(vertex_z, vertex_r);
        }
        
        n_brems++;
      }
      to_print += n_brems > 0 ? Form(" and has emitted %u bremsstrahlung photons\n", n_brems) : " and has not emitted any bremsstrahlung photons\n";
    }
    // std::cout << to_print << std::endl;
    if (primary_electron_pT_initial > 0) {
      float relative_total_brem_pT_loss = total_brem_pT_loss / primary_electron_pT_initial;
      brem_pT_loss_hist->Fill(relative_total_brem_pT_loss);
    }
    brem_photon_count_hist->Fill(n_brems);
    if (primary_electron_pT_initial > 0) { // avoid filling if we didn't find the primary electron for some reason
      // printf("PT of primary electron: %.3f GeV\n", primary_electron_pT_initial);
      brem_prob_wrt_pT_hist->Fill(n_brems > 0 ? 1 : 0, primary_electron_pT_initial);
      nPhotons_wrt_pT_profile->Fill(primary_electron_pT_initial, n_brems);
    }
  }

  // Normalize histograms by total number of events
  unsigned long n_events = particles_tree->GetEntries();
  if (n_events > 0) {
    total_pT_loss_hist->Scale(1.0 / n_events);
    brem_pT_loss_hist->Scale(1.0 / n_events);
    brem_pT_loss_hist_gen0->Scale(1.0 / n_events);
    brem_photon_count_hist->Scale(1.0 / n_events);
  }

  TCanvas* canvas = new TCanvas("canvas", "Bremsstrahlung Studies", 1200, 1600);
  canvas->Divide(2,4);

  // all axes are log scaled for visibility

  // plot both total pT loss histograms
  canvas->cd(1);
  // total_pT_loss_hist->SetLineColor(kBlack);
  // total_pT_loss_hist->SetLineWidth(2);
  brem_pT_loss_hist->SetLineColor(kRed);
  brem_pT_loss_hist->SetLineWidth(2);
  brem_pT_loss_hist_gen0->SetLineColor(kBlue);
  brem_pT_loss_hist_gen0->SetLineWidth(2);
  // total_pT_loss_hist->Draw("HIST");
  brem_pT_loss_hist->Draw("HIST");
  brem_pT_loss_hist_gen0->Draw("HIST SAME");

  // put legend in top left corner & remove edges
  TLegend* leg_total = new TLegend(0.15, 0.65, 0.4, 0.85);
  leg_total->SetBorderSize(0);
  // leg_total->AddEntry(total_pT_loss_hist, "Total pT loss", "l");
  leg_total->AddEntry(brem_pT_loss_hist, "Total Bremsstrahlung pT loss", "l");
  leg_total->AddEntry(brem_pT_loss_hist_gen0, "First Bremsstrahlung pT loss", "l");
  leg_total->Draw();

  gPad->SetLogy(1);
  gPad->SetLogx(1);

  canvas->cd(2);
  // plot average number of photons wrt pT profile
  gPad->SetLogx(1);
  nPhotons_wrt_pT_profile->Draw("PE");

  // plot bremsstrahlung photon count histogram
  canvas->cd(3);
  brem_photon_count_hist->SetLineColor(kBlack);
  brem_photon_count_hist->SetLineWidth(2);
  brem_photon_count_hist->Draw("HIST");

  gPad->SetLogy(1);
  // gPad->SetLogx(1);

  // plot bremsstrahlung probability wrt pT histogram (with errors)
  canvas->cd(4);
  brem_prob_wrt_pT_hist->Draw("AP");
  gPad->Update();  // Force ROOT to paint and create the graph
  // y range from 0 to 1 for probability TEfficency histograms
  gPad->SetLogx(1);
  brem_prob_wrt_pT_hist->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);

  // plot z-r position of first generation bremsstrahlung photon vertices
  canvas->cd(5);
  gPad->SetLogz(1);
  gen1_photon_vertex_rz_hist->Draw("COLZ");

  // plot zoomed-in z-r position of first generation bremsstrahlung photon vertices
  canvas->cd(6);
  gPad->SetLogz(1);
  gen1_photon_vertex_rz_hist_zoom->Draw("COLZ");

  // plot z-r position for high pT loss (>= 50%)
  canvas->cd(7);
  gPad->SetLogz(1);
  gen1_photon_vertex_rz_hist_highPtLoss->Draw("COLZ");

  // plot z-r position for medium pT loss (10-50%)
  canvas->cd(8);
  gPad->SetLogz(1);
  gen1_photon_vertex_rz_hist_medPtLoss->Draw("COLZ");

  canvas->SaveAs(Form("pgun_studies/figures/bremsstrahlung/%s.pdf", input_dir.Data()));
}