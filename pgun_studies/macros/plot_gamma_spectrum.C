#include <iostream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

void plot_gamma_spectrum(const TString input_data_type = "gridtriplet_2T_geo_oct25_geant4_digiCut_cotTheta2Cut",
                         const float max_eta_abs=1.0) {
  gStyle->SetOptStat(0);  // disable stats box
  TString path_prefix = input_data_type;
  TString output_prefix = TString(TString(input_data_type) + Form("_eta_%.1f", max_eta_abs));

  const std::vector<TString> pTs = {"100M", "500M", "1G", "5G", "10G"};
  const std::vector<float> pT_values_gev = {0.1, 0.5, 1.0, 5.0, 10.0};
  const std::vector<TString> pT_labels = {"100 MeV/c", "500 MeV/c", "1 GeV/c", "5 GeV/c", "10 GeV/c"};
  const std::vector<Int_t> colours = {kRed, kBlue, kGreen+3, kMagenta, kBlack};

  std::vector<TH1D*> gamma_hists_pT;
  std::vector<TH1D*> gamma_hists_p;
  std::vector<TH1D*> gen1_gamma_hists_pT;
  std::vector<TH1D*> gen1_gamma_hists_p;
  std::vector<TEfficiency*> eff_hists_wrt_gamma_p;
  std::vector<TEfficiency*> eff_hists_wrt_gamma_p_gen1;

  for (unsigned i_pT = 0; i_pT < pTs.size(); i_pT++) {
    const TString pT_string = pTs[i_pT];
    TString dir_name = Form("data/%s/electron/%s", input_data_type.Data(), pT_string.Data());
    TString particles_filename = dir_name + "/particles_simulation.root";
    
    TFile* particles_file = TFile::Open(particles_filename);
    if (!particles_file || particles_file->IsZombie()) {
      std::cout << "Warning: Could not open " << particles_filename << std::endl;
      continue;
    }
    TFile* perf_ambi_file = TFile::Open(dir_name + "/performance_finding_ambi.root");
    if (!perf_ambi_file || perf_ambi_file->IsZombie()) {
      std::cout << "Warning: Could not open " << dir_name + "/performance_finding_ambi.root" << std::endl;
      continue;
    }
    
    TTree* particles_tree = (TTree*)particles_file->Get("particles");
    TTree* perf_ambi_tree = (TTree*)perf_ambi_file->Get("matchingdetails");

    std::vector<int>* pdg = nullptr;
    std::vector<float>* eta = nullptr;
    std::vector<float>* pT = nullptr;
    std::vector<float>* p = nullptr;
    std::vector<unsigned>* generation = nullptr;

    std::vector<bool>* wasMatched = nullptr;
    std::vector<unsigned>* generation_matched = nullptr;


    particles_tree->SetBranchAddress("particle_type", &pdg);
    particles_tree->SetBranchAddress("eta", &eta);
    particles_tree->SetBranchAddress("pt", &pT);
    particles_tree->SetBranchAddress("p", &p);
    particles_tree->SetBranchAddress("generation", &generation);

    perf_ambi_tree->SetBranchAddress("matched", &wasMatched);
    perf_ambi_tree->SetBranchAddress("particle_id_generation", &generation_matched);

    // Create logarithmic binning
    const int nBins = 30;
    double logBins[nBins + 1];
    double effXBins[nBins + 1]; // for efficiency histograms (p fraction)
    double xMax = pT_values_gev[i_pT];
    double xMin = xMax / 10000.0;
    double logMin = TMath::Log10(xMin);
    double logMax = TMath::Log10(xMax);
    double logStep = (logMax - logMin) / nBins;
    for (int i = 0; i <= nBins; i++) {
      logBins[i] = TMath::Power(10, logMin + i * logStep);
      effXBins[i] = logBins[i] / xMax;
    }

    TH1D* gamma_hist_pT = new TH1D(Form("gamma_hist_%s", pT_string.Data()),
                                   Form("Gamma Spectrum for %s electrons (|#eta| < %.1f); pT (GeV/c); Counts",
                                         pT_string.Data(), max_eta_abs),
                                         nBins, logBins);
    TH1D* gen1_gamma_hist_pT = new TH1D(Form("gen1_gamma_hist_%s", pT_string.Data()),
                                        Form("Generation 1 Gamma Spectrum for %s electrons (|#eta| < %.1f); pT (GeV/c); Counts",
                                             pT_string.Data(), max_eta_abs),
                                        nBins, logBins);
    TH1D* gamma_hist_p = new TH1D(Form("gamma_hist_p_%s", pT_string.Data()),
                                  Form("Gamma Spectrum for %s electrons (|#eta| < %.1f); p (GeV/c); Counts",
                                        pT_string.Data(), max_eta_abs),
                                  nBins, logBins);
    TH1D* gen1_gamma_hist_p = new TH1D(Form("gen1_gamma_hist_p_%s", pT_string.Data()),
                                       Form("Generation 1 Gamma Spectrum for %s electrons (|#eta| < %.1f); p (GeV/c); Counts",
                                             pT_string.Data(), max_eta_abs),
                                       nBins, logBins);

    TEfficiency* eff_hist_wrt_gamma_p =
      new TEfficiency(Form("eff_hist_%s", pT_string.Data()),
                      Form("Efficiency of %s electrons wrt total #gamma p; p fraction; Efficiency",
                      pT_string.Data()), nBins, effXBins);
    TEfficiency* eff_hist_wrt_gamma_p_gen1 =
      new TEfficiency(Form("eff_hist_gen1_%s", pT_string.Data()),
                      Form("Efficiency of %s electrons wrt generation 1 #gamma p; p fraction; Efficiency",
                      pT_string.Data()), nBins, effXBins);

    for (unsigned i_ev = 0; i_ev < particles_tree->GetEntries(); i_ev++) {
      particles_tree->GetEntry(i_ev);
      perf_ambi_tree->GetEntry(i_ev);

      // get i_part_perf to point to the original electron
      unsigned i_part_perf = 0;
      while (i_part_perf < generation_matched->size()) {
        if (generation_matched->at(i_part_perf) == 0)
          break;
        i_part_perf++;
      }
      // check if electron was matched: make sure it is trackable
      // and hence written at all (to avoid segfaulting)
      bool electron_matched = false;
      if (wasMatched->size() > i_part_perf) {
        electron_matched = wasMatched->at(i_part_perf);
      }

      for (unsigned i_part = 0; i_part < pdg->size(); i_part++) {
        if (std::abs(eta->at(i_part)) > max_eta_abs)
          continue;  // skip particles outside eta range
        if (pdg->at(i_part) == 22) {  // photon
          gamma_hist_pT->Fill(pT->at(i_part));
          gamma_hist_p->Fill(p->at(i_part));
          if (generation->at(i_part) == 1) {
            gen1_gamma_hist_pT->Fill(pT->at(i_part));
            gen1_gamma_hist_p->Fill(p->at(i_part));
            eff_hist_wrt_gamma_p_gen1->Fill(electron_matched,
                                            p->at(i_part) / pT_values_gev[i_pT]);
          }
        }
      }  // loop over particles
    }  // loop over events
    gamma_hists_pT.push_back(gamma_hist_pT);
    gamma_hists_p.push_back(gamma_hist_p);
    gen1_gamma_hists_pT.push_back(gen1_gamma_hist_pT);
    gen1_gamma_hists_p.push_back(gen1_gamma_hist_p);
    eff_hists_wrt_gamma_p.push_back(eff_hist_wrt_gamma_p);
    eff_hists_wrt_gamma_p_gen1.push_back(eff_hist_wrt_gamma_p_gen1);
  }  // loop over pTs

  // Now plot gamma spectra
  TCanvas* canvas = new TCanvas("canvas_gamma_spectra", "Gamma Spectra", 1500, 1300);
  canvas->Divide(2,2);
  // pT spectra for all gammas
  canvas->cd(1);
  unsigned y_max;  // for setting y axis range
  for (unsigned i_pT = 0; i_pT < gamma_hists_pT.size(); i_pT++) {
    if (gamma_hists_pT[i_pT]->GetMaximum() > y_max) {
      y_max = gamma_hists_pT[i_pT]->GetMaximum();
    }
  }
  y_max = static_cast<unsigned>(y_max * 1.1);  // 10% margin

  for (unsigned i_pT = 0; i_pT < gamma_hists_pT.size(); i_pT++) {
    unsigned i_pT_to_plot = gamma_hists_pT.size() - 1 - i_pT;  // reverse order for better visibility
    gamma_hists_pT[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
    gamma_hists_pT[i_pT_to_plot]->SetLineWidth(2);
    if (i_pT == 0) {
      gamma_hists_pT[i_pT_to_plot]->Draw("HIST");
      gamma_hists_pT[i_pT_to_plot]->SetTitle(Form("#gamma Spectrum of e^{-} brem at various e^{-} p_{T} (|#eta| < %.1f); p_{T} (GeV/c); Counts", max_eta_abs));
      gamma_hists_pT[i_pT_to_plot]->GetYaxis()->SetRangeUser(1, y_max);
    } else {
      gamma_hists_pT[i_pT_to_plot]->Draw("HIST SAME");
    }
  }
  // same legend for all four plots in top right
  TLegend* leg = new TLegend(0.55, 0.65, 0.8, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (unsigned i_pT = 0; i_pT < pTs.size(); i_pT++) {
    leg->AddEntry(gamma_hists_pT[i_pT], TString("e^{-} p_{T} = ") + pT_labels[i_pT], "l");
  }
  leg->Draw();
  // gPad->SetLogy(1);
  gPad->SetLogx(1);

  // p spectra for all gammas
  canvas->cd(2);
  y_max = 0;  // reset for setting y axis range
  for (unsigned i_pT = 0; i_pT < gamma_hists_p.size(); i_pT++) {
    if (gamma_hists_p[i_pT]->GetMaximum() > y_max) {
      y_max = gamma_hists_p[i_pT]->GetMaximum();
    }
  }
  y_max = static_cast<unsigned>(y_max * 1.1);

  for (unsigned i_pT = 0; i_pT < gamma_hists_p.size(); i_pT++) {
    unsigned i_pT_to_plot = gamma_hists_p.size() - 1 - i_pT;  // reverse order for better visibility
    gamma_hists_p[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
    gamma_hists_p[i_pT_to_plot]->SetLineWidth(2);
    if (i_pT == 0) {
      gamma_hists_p[i_pT_to_plot]->Draw("HIST");
      gamma_hists_p[i_pT_to_plot]->SetTitle(Form("#gamma Spectrum of e^{-} brem at various e^{-} p_{T} (|#eta| < %.1f); p (GeV/c); Counts", max_eta_abs));
      gamma_hists_p[i_pT_to_plot]->GetYaxis()->SetRangeUser(1, y_max);
    } else {
      gamma_hists_p[i_pT_to_plot]->Draw("HIST SAME");
    }
  }
  leg->Draw();
  // gPad->SetLogy(1);
  gPad->SetLogx(1);

  // pT spectra for gen 0 gammas
  y_max = 0;  // reset for setting y axis range
  for (unsigned i_pT = 0; i_pT < gen1_gamma_hists_pT.size(); i_pT++) {
    if (gen1_gamma_hists_pT[i_pT]->GetMaximum() > y_max) {
      y_max = gen1_gamma_hists_pT[i_pT]->GetMaximum();
    }
  }
  y_max = static_cast<unsigned>(y_max * 1.1);

  canvas->cd(3);
  for (unsigned i_pT = 0; i_pT < gen1_gamma_hists_pT.size(); i_pT++) {
    unsigned i_pT_to_plot = gen1_gamma_hists_pT.size() - 1 - i_pT;  // reverse order for better visibility
    gen1_gamma_hists_pT[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
    gen1_gamma_hists_pT[i_pT_to_plot]->SetLineWidth(2);
    if (i_pT == 0) {
      gen1_gamma_hists_pT[i_pT_to_plot]->Draw("HIST");
      gen1_gamma_hists_pT[i_pT_to_plot]->SetTitle(Form("Gen 1 #gamma Spectrum of e^{-} brem at various e^{-} p_{T} (|#eta| < %.1f); p_{T} (GeV/c); Counts", max_eta_abs));
      gen1_gamma_hists_pT[i_pT_to_plot]->GetYaxis()->SetRangeUser(1, y_max);
    } else {
      gen1_gamma_hists_pT[i_pT_to_plot]->Draw("HIST SAME");
    }
  }
  leg->Draw();
  // gPad->SetLogy(1);
  gPad->SetLogx(1);

  // p spectra for gen 0 gammas
  canvas->cd(4);
  y_max = 0;  // reset for setting y axis range
  for (unsigned i_pT = 0; i_pT < gen1_gamma_hists_p.size(); i_pT++) {
    if (gen1_gamma_hists_p[i_pT]->GetMaximum() > y_max) {
      y_max = gen1_gamma_hists_p[i_pT]->GetMaximum();
    }
  }
  y_max = static_cast<unsigned>(y_max * 1.1);

  for (unsigned i_pT = 0; i_pT < gen1_gamma_hists_p.size(); i_pT++) {
    unsigned i_pT_to_plot = gen1_gamma_hists_p.size() - 1 - i_pT;  // reverse order for better visibility
    gen1_gamma_hists_p[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
    gen1_gamma_hists_p[i_pT_to_plot]->SetLineWidth(2);
    if (i_pT == 0) {
      gen1_gamma_hists_p[i_pT_to_plot]->Draw("HIST");
      gen1_gamma_hists_p[i_pT_to_plot]->SetTitle(Form("Gen 1 #gamma Spectrum of e^{-} brem at various e^{-} p_{T} (|#eta| < %.1f); p (GeV/c); Counts", max_eta_abs));
      gen1_gamma_hists_p[i_pT_to_plot]->GetYaxis()->SetRangeUser(1, y_max);
    } else {
      gen1_gamma_hists_p[i_pT_to_plot]->Draw("HIST SAME");
    }
  }
  leg->Draw();
  // gPad->SetLogy(1);
  gPad->SetLogx(1);

  canvas->SaveAs(Form("figures/gamma_spectra/%s.pdf", output_prefix.Data()));
  delete canvas;

  // now plot efficiency wrt gamma p
  canvas = new TCanvas("canvas_eff_gamma_p", "Efficiency wrt gamma p", 1500, 600);
  canvas->Divide(2,1);
  // efficiency wrt all gammas
  canvas->cd(1);
  // for (unsigned i_pT = 0; i_pT < eff_hists_wrt_gamma_p.size(); i_pT++) {
  //   unsigned i_pT_to_plot = eff_hists_wrt_gamma_p.size() - 1 - i_pT;  // reverse order for better visibility
  //   eff_hists_wrt_gamma_p[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
  //   eff_hists_wrt_gamma_p[i_pT_to_plot]->SetLineWidth(2);
  //   if (i_pT == 0) {
  //     eff_hists_wrt_gamma_p[i_pT_to_plot]->Draw("AP");
  //     eff_hists_wrt_gamma_p[i_pT_to_plot]->SetTitle(
  //       Form("Electron Efficiency wrt all #gamma p at various e^{-} p_{T} (|#eta| < %.1f); p fraction; Efficiency", max_eta_abs));
  //     eff_hists_wrt_gamma_p[i_pT_to_plot]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  //   } else {
  //     eff_hists_wrt_gamma_p[i_pT_to_plot]->Draw("AP SAME");
  //   }
  // }
  // leg->Draw();
  // gPad->SetLogx(1);

  canvas->cd(2);
  // efficiency wrt gen 0 gammas
  for (unsigned i_pT = 0; i_pT < eff_hists_wrt_gamma_p_gen1.size(); i_pT++) {
    unsigned i_pT_to_plot = eff_hists_wrt_gamma_p_gen1.size() - 1 - i_pT;  // reverse order for better visibility
    eff_hists_wrt_gamma_p_gen1[i_pT_to_plot]->SetLineColor(colours[i_pT_to_plot]);
    eff_hists_wrt_gamma_p_gen1[i_pT_to_plot]->SetLineWidth(2);
    if (i_pT == 0) {
      eff_hists_wrt_gamma_p_gen1[i_pT_to_plot]->Draw("P");
      eff_hists_wrt_gamma_p_gen1[i_pT_to_plot]->SetTitle(
        Form("Electron Efficiency wrt gen 1 #gamma p at various e^{-} p_{T} (|#eta| < %.1f); p fraction; Efficiency", max_eta_abs));
    } else {
      eff_hists_wrt_gamma_p_gen1[i_pT_to_plot]->Draw("P SAME");
    }
  }
  gPad->Update();  // Need to update before accessing painted graph
  eff_hists_wrt_gamma_p_gen1[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  leg->Draw();
  gPad->SetLogx(1);

  canvas->SaveAs(Form("figures/gamma_spectra/%s_eff_wrt_gamma_p.pdf", output_prefix.Data()));
  delete canvas;
}