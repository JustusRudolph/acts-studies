#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>


void plot_efficiency_analysis(bool use_anti = false, const std::string input_data_type="truth_seeds") {
  
  // Set ROOT style
  gStyle->SetOptStat(0);
  
  // Define particle types and pTs
  std::vector<TString> reg_particles = {"electron", "muon", "piplus", "proton", "kaonplus"};
  std::vector<TString> reg_particles_labels = {"e^{-}", "#mu^{-}", "#pi^{+}", "p^{+}", "K^{+}"};
  std::vector<TString> anti_particles = {"positron", "antimuon", "piminus", "antiproton", "kaonminus"};
  std::vector<TString> anti_particles_labels = {"e^{+}", "#mu^{+}", "#pi^{-}", "p^{-}", "K^{-}"};
  std::vector<int> pTs_mev = {50, 60, 70, 80, 90, 100, 200, 500, 1000, 10000, 100000};
  // to have central bin values correct for plotting
  std::array<double, 12> pT_bin_edges_mev = {45., 55., 65., 75., 85., 95., 105., 295., 705., 1295., 18705., 181295.};
  std::vector<TString> pT_strings = {"50M", "60M", "70M", "80M", "90M", "100M", "200M", "500M", "1G", "10G", "100G"};
  std::vector<TString> pT_labels = {"50 MeV/c", "60 MeV/c", "70 MeV/c", "80 MeV/c", "90 MeV/c", "100 MeV/c",
                                    "200 MeV/c", "500 MeV/c", "1 GeV/c", "10 GeV/c", "100 GeV/c"};
  std::vector<unsigned int> energies_to_plot = {0, 5, 8, 9}; // indices of pT_strings

  std::vector<TString> particles;
  std::vector<TString> particle_labels;
  if (use_anti) {
    particles = anti_particles;
    particle_labels = anti_particles_labels;
  } else {
    particles = reg_particles;
    particle_labels = reg_particles_labels;
  }

  // Colours for different pT_strings
  std::vector<Int_t> colours_per_particle = {kRed, kBlue, kGreen+2, kYellow+1, kMagenta};
  
  // Store histograms for each particle and energy
  std::map<TString, std::map<TString, TEfficiency*>> efficiencyHistosEta;
  std::map<TString, std::map<TString, TEfficiency*>> efficiencyHistosPhi;
  std::map<TString, std::map<TString, float>> efficiencies;
  std::vector<TH1D*> eff_hists(particles.size());
  // std::map<TString, std::map<TString, float>> fake_rates;

  // First, try to load all available files and extract efficiency histograms
  for (size_t i = 0; i < particles.size(); i++) {
    const auto& particle = particles[i];
    eff_hists[i] = new TH1D(Form("eff_hist_%s", particle.Data()), Form("Efficiency vs p_{T} for %s", particle.Data()),
                            pT_bin_edges_mev.size() - 1, pT_bin_edges_mev.data());
    for (size_t j = 0; j < pT_strings.size(); j++) {
      const auto& energy = pT_strings[j];
      TString filename = Form("data/%s/%s/%s/performance_finding_ambi.root", 
                             input_data_type.c_str(), particle.Data(), energy.Data());
      
      TFile* file = TFile::Open(filename);
      if (!file || file->IsZombie()) {
        std::cout << "Warning: Could not open " << filename << std::endl;
        continue;
      }
      
      TVectorT<float>* eff_particles = (TVectorT<float>*)file->Get("eff_particles");
      efficiencies[particle][energy] = (*eff_particles)[0];
      eff_hists[i]->SetBinContent(j+1, (*eff_particles)[0]);
      // printf("WRITING GRAPH: Efficiency for %s at %s: %f\n", particle.Data(), energy.Data(), (*eff_particles)[0]);
      
      // Look for efficiency histograms in the file
      TEfficiency* hist_eta = (TEfficiency*)file->Get("trackeff_vs_eta");
      TEfficiency* hist_eta_clone = (TEfficiency*)hist_eta->Clone();
      TEfficiency* hist_phi = (TEfficiency*)file->Get("trackeff_vs_phi");
      TEfficiency* hist_phi_clone = (TEfficiency*)hist_phi->Clone();
      hist_eta_clone->SetName(Form("eff_vs_eta_%s_%s", particle.Data(), energy.Data()));
      hist_phi_clone->SetName(Form("eff_vs_phi_%s_%s", particle.Data(), energy.Data()));

      if (!hist_eta) {
        std::cout << "Warning: Could not find efficiency histogram in " << filename << std::endl;
        continue;
      }
      
      // Clone the histogram and set properties
      efficiencyHistosEta[particle][energy] = hist_eta_clone;
      efficiencyHistosEta[particle][energy]->SetTitle(Form("%s %s Efficiency vs #eta", particle.Data(), energy.Data()));
      efficiencyHistosEta[particle][energy]->SetLineColor(colours_per_particle[i]);
      efficiencyHistosEta[particle][energy]->SetMarkerColor(colours_per_particle[i]);

      efficiencyHistosPhi[particle][energy] = hist_phi_clone;
      efficiencyHistosPhi[particle][energy]->SetTitle(Form("%s %s Efficiency vs #phi", particle.Data(), energy.Data()));
      efficiencyHistosPhi[particle][energy]->SetLineColor(colours_per_particle[i]);
      efficiencyHistosPhi[particle][energy]->SetMarkerColor(colours_per_particle[i]);
      
      file->Close();
      delete file;
    }
  }
  
  // Create canvas with 4 rows, 2 cols
  TCanvas* canvas = new TCanvas("canvas", "Efficiency Analysis", 1200, 1600);
  canvas->Divide(2, 4);
  
  // Plot 1: Efficiency wrt pT for all particles
  canvas->cd(1);
  TLegend* leg1 = new TLegend(0.45, 0.15, 0.75, 0.35);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  eff_hists[0]->SetTitle("Efficiency vs p_{T} for various particles; p_{T} (MeV/c); Efficiency");
  eff_hists[0]->GetXaxis()->SetTitleOffset(1.2);  // avoid axis title overwriting axis

  for (size_t i = 0; i < particles.size(); i++) {
    eff_hists[i]->SetMarkerColor(colours_per_particle[i]);
    eff_hists[i]->SetLineColor(colours_per_particle[i]);
    eff_hists[i]->SetMarkerStyle(2);
    eff_hists[i]->SetMarkerSize(1);
    eff_hists[i]->SetLineStyle(2);
    eff_hists[i]->SetLineWidth(1);
    leg1->AddEntry(eff_hists[i], particle_labels[i], "lp");
    if (i == 0) {
      eff_hists[i]->Draw("LP");
      eff_hists[i]->GetYaxis()->SetRangeUser(0., 1.05);
      eff_hists[i]->GetXaxis()->SetRangeUser(50, 5000);
    } else {
      eff_hists[i]->Draw("LP SAME");
    }
  }
  gPad->SetLogx(1);
  gPad->Update();

  leg1->Draw();

  // Energy comparisons (electron, muon, pion)
  TLegend* leg2 = new TLegend(0.45, 0.15, 0.75, 0.35);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(efficiencyHistosEta[particles[0]][pT_strings[0]], particle_labels[0], "lp");
  leg2->AddEntry(efficiencyHistosEta[particles[1]][pT_strings[0]], particle_labels[1], "lp");
  leg2->AddEntry(efficiencyHistosEta[particles[2]][pT_strings[0]], particle_labels[2], "lp");
  leg2->AddEntry(efficiencyHistosEta[particles[3]][pT_strings[0]], particle_labels[3], "lp");
  leg2->AddEntry(efficiencyHistosEta[particles[4]][pT_strings[0]], particle_labels[4], "lp");

  for (size_t i = 0; i < energies_to_plot.size(); i++) {
    canvas->cd(i+2);
    // printf("Plotting %s\n", pT_labels[i].Data());
    unsigned int energy_index = energies_to_plot[i];

    efficiencyHistosEta[particles[0]][pT_strings[energy_index]]->SetTitle(
      Form("%s Efficiency Comparison; #eta; Efficiency", pT_labels[energy_index].Data()));
    efficiencyHistosEta[particles[0]][pT_strings[energy_index]]->Draw("AP");
    efficiencyHistosEta[particles[1]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosEta[particles[2]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosEta[particles[3]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosEta[particles[4]][pT_strings[energy_index]]->Draw("P SAME");
    gPad->Update();
    
    if (pTs_mev[energy_index] < 1000) {
      efficiencyHistosEta[particles[0]][pT_strings[energy_index]]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);
    } else {
      efficiencyHistosEta[particles[0]][pT_strings[energy_index]]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.5, 1.05);
    }

    leg2->Draw();
    canvas->Update();
  }
  canvas->Update();
  // Save the plot
  TString output_pdf = Form("figures/efficiency_analysis/%s_%s.pdf", use_anti ? "anti" : "reg", input_data_type.c_str());
  canvas->SaveAs(output_pdf);

  TCanvas* canvas_phi = new TCanvas("canvas_phi", "Efficiency Analysis Phi", 1600, 600 * (energies_to_plot.size() / 2));
  canvas_phi->Divide(2, energies_to_plot.size() / 2);
  for (size_t i = 0; i < energies_to_plot.size(); i++) {
    canvas_phi->cd(i+1);
    // printf("Plotting %s\n", pT_labels[i].Data());
    unsigned int energy_index = energies_to_plot[i];

    efficiencyHistosPhi[particles[0]][pT_strings[energy_index]]->SetTitle(
      Form("%s Efficiency Comparison; #phi; Efficiency", pT_labels[energy_index].Data()));
    efficiencyHistosPhi[particles[0]][pT_strings[energy_index]]->Draw("AP");
    efficiencyHistosPhi[particles[1]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosPhi[particles[2]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosPhi[particles[3]][pT_strings[energy_index]]->Draw("P SAME");
    efficiencyHistosPhi[particles[4]][pT_strings[energy_index]]->Draw("P SAME");
    gPad->Update();
    
    if (pTs_mev[energy_index] < 1000) {
      efficiencyHistosPhi[particles[0]][pT_strings[energy_index]]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);
    } else {
      efficiencyHistosPhi[particles[0]][pT_strings[energy_index]]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);
    }

    leg2->Draw();
  }
  TString output_pdf_phi = Form("figures/efficiency_analysis/phi_%s_%s.pdf", use_anti ? "anti" : "reg", input_data_type.c_str());
  canvas_phi->SaveAs(output_pdf_phi);

  return;
} 