{
  const float zVtx = 0.0;
  
  TCanvas *c1 = new TCanvas("c1","c1: geometry", 800, 600);
  auto h1 = gPad->DrawFrame(-2.5,0,2.8,1);
  
  h1->GetXaxis()->SetTitle("x (m)");
  h1->GetYaxis()->SetTitle("r (m)");

  // geometry in mm

  //const float MLBz = 1288; // default
  const float MLBz = 12*128.8; // 12 modules
  const int nMLB = 5;
  const float MLBr[nMLB] = {70, 90, 120, 200, 300}; // last layer is formally OT project
  const float MLDrmin = 100;
  const float MLDrmax = 350;
  const int nMLD = 3;
  //const float MLDz[nMLD] = {770, 1000, 1220}; 
  //const float MLDz[nMLD] = {860, 1000, 1220}; // barrel 11 modules + 15 cm gap
  const float MLDz[nMLD] = {920, 1000, 1220}; // barrel 11 modules + 15 cm gap
  
  const float OTBz = 2*1288; // default
  //const float OTBz = 22*128.8;
  const int nOTB = 3;
  const float OTBr[nOTB] = {450, 600, 800};
  const float OTDrmin = 200;
  const float OTDrmax = 680;
  const int nOTD = 3;
  //const float OTDz[nOTD] = {1500, 1800, 2200}; // default
  const float OTDz[nOTD] = {1440, 1800, 2200}; // gap = 15 cm
  //const float OTDz[nOTD] = {1570, 1800, 2200}; // 11 modules gap = 15 cm
  
  const Float_t kScale = 1e-3; // mm to m for visualisation

  TLine *line = new TLine();
  line->SetLineWidth(3);
  for (Int_t iL = 0; iL < nMLB; iL++) {
    line->DrawLine(-MLBz/2*kScale, MLBr[iL]*kScale, MLBz/2*kScale, MLBr[iL]*kScale);
  }
  for (Int_t iD = 0; iD < nMLD; iD++) {
    line->DrawLine(-MLDz[iD]*kScale, MLDrmin*kScale, -MLDz[iD]*kScale, MLDrmax*kScale);
    line->DrawLine(MLDz[iD]*kScale, MLDrmin*kScale, MLDz[iD]*kScale, MLDrmax*kScale);
  }
 
  for (Int_t iL = 0; iL < nOTB; iL++) {
    line->DrawLine(-OTBz/2*kScale, OTBr[iL]*kScale, OTBz/2*kScale, OTBr[iL]*kScale);
  }
  for (Int_t iD = 0; iD < nOTD; iD++) {
    line->DrawLine(-OTDz[iD]*kScale, OTDrmin*kScale, -OTDz[iD]*kScale, OTDrmax*kScale);
    line->DrawLine(OTDz[iD]*kScale, OTDrmin*kScale, OTDz[iD]*kScale, OTDrmax*kScale);
  }

  float rapVals[] = {0, 1.3, 1.5, 1.7, 1.9, 2.5};
  // float rapVals[] = {0, 1.3, 1.68, 2.5};
  int nRap = sizeof(rapVals)/sizeof(rapVals[0]);
  float rrap = 0.9;
  float zrap = 2.4;
  line->SetLineWidth(1);
  line->SetLineColor(kGray+2);
  TLatex *ltx = new TLatex();
  ltx->SetTextFont(42);
  ltx->SetTextSize(0.025);
  
  for (int iRap = 0; iRap < nRap; iRap++) {
    float rEnd = rrap;
    float zEnd = zVtx + rEnd/TMath::Tan(2*TMath::ATan(TMath::Exp(-rapVals[iRap])));
    const char* label = Form("#eta=%.2f",rapVals[iRap]);
    if (zEnd > zrap) {
      zEnd = zrap;
      rEnd = (zEnd - zVtx) * TMath::Tan(2*TMath::ATan(TMath::Exp(-rapVals[iRap])));
      ltx->SetTextAlign(12);
      ltx->DrawLatex(zEnd+0.02,rEnd,label);
    }
    else {
      ltx->SetTextAlign(21);
      ltx->DrawLatex(zEnd,rEnd+0.02,label);
    }
    line->DrawLine(zVtx,0,zEnd,rEnd);
  }
  // save the canvas
  c1->SaveAs("figures/layout/geo_layout_eta_lines.pdf");
}
