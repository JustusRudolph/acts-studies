#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void geo_efficiency_plots(
  const std::vector<TString> pathBases = {"geo_staves_pi_1GeV_eta14-20"},
  const bool onSTBC=false,
  const unsigned nMinHits=3,
  const float etaMin=-2.5,
  const float etaMax=2.5
)
{
  const TString base = onSTBC ? "/data/alice/jrudolph/" : "/home/justus/projects/";
  const TString actso2_output_base    = base + "alice/ACTSO2/output/";
  const TString geo_studies_base = base + "alice/acts-studies/geo_studies/";
  const TString rootOutFile = geo_studies_base + "histos/efficiencies/geo_efficiency_plots.root";
  const TString pdfOutFile = geo_studies_base + "figures/efficiencies/";

  // Plot 1: mean number of hits vs eta (100 bins, etaMin–etaMax)
  TProfile* pHitsAll   = new TProfile("pHitsAll",   ";;Mean number of hits", 100, etaMin, etaMax);
  TProfile* pHitsReco  = new TProfile("pHitsReco",  ";;Mean number of hits", 100, etaMin, etaMax);
  TProfile* pHitsNReco = new TProfile("pHitsNReco", ";;Mean number of hits", 100, etaMin, etaMax);

  // Plot 2: 1D efficiency vs eta (100 bins, etaMin–etaMax)
  auto* effEta = new TEfficiency("effEta", ";#eta;Efficiency",
                                 100, etaMin, etaMax);

  // Plot 3: 2D efficiency eta-phi (etamMin-etaMax, full phi)
  auto* eff2D = new TEfficiency("eff2D", ";#eta;#phi;Efficiency",
                                100, etaMin, etaMax, 100, -TMath::Pi(), TMath::Pi());

  // Plot 4: 2D efficiency eta-pT (100 bins eta, etaMin–etaMax)
  auto* effEtaPt = new TEfficiency("effEtaPt", ";#eta;p_{T} (GeV/c);Efficiency",
                                   100, etaMin, etaMax, 100, 0.5, 5.0);

  // --- Event loop over all input files ---
  unsigned nFillsTotal{0};
  for (const TString& pathBase : pathBases) {
    unsigned nFills{0};
    std::cout << "Processing " << pathBase << "..." << std::endl;

    const TString inFile = actso2_output_base + pathBase + "/particles_matched.root";
    TFile* fIn = TFile::Open(inFile);
    if (!fIn || fIn->IsZombie()) { std::cerr << "Cannot open " << inFile << " — skipping" << std::endl; continue; }

    TTree* t = (TTree*)fIn->Get("particles");
    if (!t) { std::cerr << "particles tree not found in " << inFile << " — skipping" << std::endl; fIn->Close(); continue; }

    std::vector<float>* eta = nullptr;
    std::vector<float>* phi = nullptr;
    std::vector<float>* pt = nullptr;
    std::vector<int>* nhits = nullptr;
    std::vector<unsigned>* generation = nullptr;
    std::vector<std::vector<unsigned int>>* matchedIdxs = nullptr;

    t->SetBranchAddress("eta",                &eta);
    t->SetBranchAddress("phi",                &phi);
    t->SetBranchAddress("pt",                 &pt);
    t->SetBranchAddress("number_of_hits",     &nhits);
    t->SetBranchAddress("generation",         &generation);
    t->SetBranchAddress("matched_track_idxs", &matchedIdxs);

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
      t->GetEntry(i);

      for (unsigned i_part = 0; i_part < eta->size(); i_part++) {
        if (generation->at(i_part) != 0 ||
            nhits->at(i_part) < nMinHits) continue;

        const bool wasReco = !matchedIdxs->at(i_part).empty();
        nFills++;
        nFillsTotal++;

        pHitsAll->Fill(eta->at(i_part), nhits->at(i_part));
        if (wasReco) pHitsReco->Fill(eta->at(i_part), nhits->at(i_part));
        else         pHitsNReco->Fill(eta->at(i_part), nhits->at(i_part));

        effEta->Fill(wasReco, eta->at(i_part));
        eff2D->Fill(wasReco, eta->at(i_part), phi->at(i_part));
        effEtaPt->Fill(wasReco, eta->at(i_part), pt->at(i_part));
      }
    }
    fIn->Close();
    std::cout << "Number of fills in " << pathBase << ": " << nFills << std::endl;
  }
  std::cout << "Total number of fills: " << nFillsTotal << std::endl;

  // --- Canvas 1: mean hits vs eta ---
  TCanvas* c1 = new TCanvas("c1", "Mean hits vs #eta", 800, 600);
  pHitsAll->SetTitle("Mean number of hits vs #eta;#eta;Mean number of hits");
  pHitsAll->SetLineColor(kBlack); pHitsAll->SetLineWidth(2);
  pHitsReco->SetLineColor(kGreen + 2); pHitsReco->SetLineWidth(2);
  pHitsNReco->SetLineColor(kRed); pHitsNReco->SetLineWidth(2);
  pHitsAll->SetStats(false);
  pHitsAll->Draw("E1");
  pHitsReco->Draw("E1 same");
  pHitsNReco->Draw("E1 same");
  pHitsAll->GetYaxis()->SetRangeUser(6, 15);
  auto* leg1 = new TLegend(0.12, 0.72, 0.38, 0.88);
  leg1->AddEntry(pHitsAll,   "All",      "ep");
  leg1->AddEntry(pHitsReco,  "Reco",     "ep");
  leg1->AddEntry(pHitsNReco, "Not reco", "ep");
  leg1->Draw();

  // --- Canvas 2: 1D efficiency vs eta ---
  TCanvas* c2 = new TCanvas("c2", "Efficiency vs #eta", 800, 600);
  effEta->Draw("AP");
  gPad->Update();
  effEta->GetPaintedGraph()->GetXaxis()->SetTitle("#eta");
  effEta->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
  effEta->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.05);
  gPad->Update();

  // --- Canvas 3: 2D efficiency eta-phi ---
  TCanvas* c3 = new TCanvas("c3", "2D Efficiency #eta-#phi", 900, 700);
  c3->SetRightMargin(0.13);
  eff2D->Draw("colz");
  gPad->Update();

  // --- Canvas 4: 2D efficiency eta-pT ---
  TCanvas* c4 = new TCanvas("c4", "2D Efficiency #eta-p_{T}", 900, 700);
  c4->SetRightMargin(0.13);
  effEtaPt->Draw("colz");
  gPad->Update();

  // --- Save ---
  c1->Print(pdfOutFile + "hits_wrt_eta_min" + TString(std::to_string(nMinHits)) + ".pdf");
  c2->Print(pdfOutFile + "eff_eta_narrow_min" + TString(std::to_string(nMinHits)) + ".pdf");
  c3->Print(pdfOutFile + "eff_eta_phi_min" + TString(std::to_string(nMinHits)) + ".pdf");
  c4->Print(pdfOutFile + "eff_eta_pt_min" + TString(std::to_string(nMinHits)) + ".pdf");
  
  TFile* fOut = TFile::Open(rootOutFile, "RECREATE");
  pHitsAll->Write(); pHitsReco->Write(); pHitsNReco->Write();
  effEta->Write("effEta");
  eff2D->Write("eff2D");
  effEtaPt->Write("effEtaPt");
  fOut->Close();

  std::cout << "Histograms saved to " << rootOutFile << std::endl;
}
