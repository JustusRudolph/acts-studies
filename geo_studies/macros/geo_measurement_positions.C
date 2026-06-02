#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>

void geo_measurement_positions(
  const TString pathToFilesFromOutput = "geo_staves_pi_1GeV_eta14-20",
  const bool useOnlyPrimaryTracks = true)
{
  const TString base           = "/home/justus/projects/alice/ACTSO2/output/";
  const TString geo_studies_base = "/home/justus/projects/alice/acts-studies/geo_studies/";
  const TString inFile         = base + pathToFilesFromOutput + "/measurements.root";
  const TString rootOutFile    = geo_studies_base + "histos/geo/geo_measurement_positions.root";
  const TString pdfOutFile     = geo_studies_base + "figures/geo/geo_measurement_positions.pdf";

  TFile* fIn = TFile::Open(inFile);
  if (!fIn || fIn->IsZombie()) { std::cerr << "Cannot open " << inFile << std::endl; return; }

  TTree* t = (TTree*)fIn->Get("measurements");
  if (!t) { std::cerr << "measurements tree not found" << std::endl; return; }

  float rec_gx, rec_gy, rec_gz;
  std::vector<unsigned int>* particles_generation = nullptr;

  t->SetBranchAddress("rec_gx",               &rec_gx);
  t->SetBranchAddress("rec_gy",               &rec_gy);
  t->SetBranchAddress("rec_gz",               &rec_gz);
  t->SetBranchAddress("particles_generation", &particles_generation);

  TH2F* hZR = new TH2F("hZR", ";z (mm);r (mm)", 200, 0, 2500, 200, 0, 1000);
  hZR->SetStats(false);

  Long64_t nEntries = t->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    t->GetEntry(i);

    if (useOnlyPrimaryTracks) {
      bool hasPrimary = false;
      for (unsigned gen : *particles_generation)
        if (gen == 0) { hasPrimary = true; break; }
      if (!hasPrimary) continue;
    }

    const float r = TMath::Sqrt(rec_gx * rec_gx + rec_gy * rec_gy);
    hZR->Fill(rec_gz, r);
  }

  TCanvas* c1 = new TCanvas("c1", "Measurement positions z vs r", 900, 700);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.08);
  c1->SetBottomMargin(0.12);
  c1->SetLogz();
  hZR->Draw("colz");
  gPad->Update();

  // --- eta reference lines (same maths as drawLayers.C) ---
  const float zMax = hZR->GetXaxis()->GetXmax();
  const float rMax = hZR->GetYaxis()->GetXmax();
  const float rapVals[] = {1.4, 1.6, 1.8, 2.0, 2.5};
  const int   nRap      = sizeof(rapVals) / sizeof(rapVals[0]);

  TLine* etaLine = new TLine();
  etaLine->SetLineWidth(1);
  etaLine->SetLineColor(kGray + 2);
  etaLine->SetLineStyle(2);

  TLatex* ltx = new TLatex();
  ltx->SetTextFont(42);
  ltx->SetTextSize(0.030);

  for (int iRap = 0; iRap < nRap; iRap++) {
    const float tanTheta = TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-rapVals[iRap])));
    // assume line clips at r = rMax first
    float rEnd = rMax;
    float zEnd = rEnd / tanTheta;
    if (zEnd > zMax) {
      // clips at z = zMax instead
      zEnd = zMax;
      rEnd = zEnd * tanTheta;
      ltx->SetTextAlign(32);
      ltx->DrawLatex(zEnd - 0.01 * zMax, rEnd + 0.02 * rMax, Form("#eta=%.1f", rapVals[iRap]));
    } else {
      ltx->SetTextAlign(23);
      ltx->DrawLatex(zEnd, rEnd - 0.045 * rMax, Form("#eta=%.1f", rapVals[iRap]));
    }
    etaLine->DrawLine(0, 0, zEnd, rEnd);
  }
  gPad->Update();

  c1->Print(pdfOutFile);

  TFile* fOut = TFile::Open(rootOutFile, "RECREATE");
  hZR->Write();
  c1->Write();
  fOut->Close();
  fIn->Close();

  std::cout << "Plot saved to " << pdfOutFile << std::endl;
}
