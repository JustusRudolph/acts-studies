#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"

#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TCanvas.h"

#include "TPad.h"
#include "TStyle.h"
#include "TChain.h"

int pmObjectCounter = 0;

void drawPM(int mStyle, float mSize, int color, bool drawLines = true)
{
  TPolyMarker3D *graph_states = (TPolyMarker3D *)gPad->GetPrimitive("TPolyMarker3D");

  if (graph_states != 0x0)
  {
    // graph_hits->SetTitle(";x, mm;y, mm");
    graph_states->SetName(Form("3D_pm_%d", pmObjectCounter++));
    graph_states->SetMarkerColor(color); // i==0 ? kRed : kBlue+1);
    graph_states->SetMarkerStyle(mStyle);
    graph_states->SetMarkerSize(mSize);
    cout << "n points in pm = " << graph_states->GetN() << endl;

    // lines, if needed
    if (drawLines)
    {
      TPolyLine3D *pl = new TPolyLine3D(0); // graph_states->GetN());
      cout << pl->GetN() << endl;
      // pl->SetName("pl_states");
      for (int i = 0; i < graph_states->GetN(); i++)
      {
        double x, y, z;
        graph_states->GetPoint(i, x, y, z);
        if (isnan(x))
            continue;
        // cout << "x=" << x << ", y=" << y << ", z=" << z << endl;
        pl->SetNextPoint(x, y, z);
      }
      pl->SetLineWidth(2);
      pl->SetLineColor(color);
      pl->Draw("same");
    }

    // leg->AddEntry(graph_states, strFileNames[i][1], "p");
  }
}

// ###########
// ###########
// ###########
void draw_hits_seeds_states_EbyE_Nov26_Re_Jan2026(TString path_to_hits,
                                                  unsigned evToDraw,
                                                  bool drawSecondaries = false)
{
  TChain *inputTreeHits;
  TChain *inputTreeStates;

  TLegend *leg = new TLegend(0.6, 0.78, 0.88, 0.89);
  leg->SetFillColor(kWhite);

  TCanvas *canv_hits = new TCanvas(Form("canv_hits_%d", evToDraw),
                                   Form("canv_hits_%d", evToDraw),
                                   1, 1, 1600, 800);
  canv_hits->Divide(2, 1);
  
  // ========== PAD 1: 3D View ==========
  canv_hits->cd(1);
  
  // setting frame + origin
  TPolyMarker3D *grFrame = new TPolyMarker3D(3);
  grFrame->SetPoint(0, -1000, -800, -800);
  grFrame->SetPoint(1, 1000, 800, 800);
  grFrame->SetName("tmpFrame");
  grFrame->SetMarkerColor(kWhite);
  TView *view = TView::CreateView(1);
  view->SetRange(-1000, -800, -800, 1000, 800, 800);
  // view->SetRange(-500, -400, -400, 500, 400, 400);
  // view->SetRange(-400, -300, -300, 400, 300, 300);
  grFrame->DrawClone("APL");

  TPolyMarker3D *grOrigin = new TPolyMarker3D(3);
  grOrigin->SetPoint(0, 0, 0, 0);
  grOrigin->SetName("grOrigin");
  grOrigin->SetMarkerColorAlpha(kGreen + 1, 0.8);
  grOrigin->SetMarkerStyle(29);
  grOrigin->SetMarkerSize(2.2);
  grOrigin->DrawClone("P");

  auto ax = new TAxis3D();
  ax->SetAxisColor(kBlack);
  ax->SetLabelColor(kBlack);
  ax->SetXTitle("z, mm");
  ax->SetYTitle("x, mm");
  ax->SetZTitle("y, mm");
  ax->SetTitleOffset(1.4);
  ax->Draw();
  gPad->Update();

  // ### hits
  inputTreeHits = new TChain("hits");
  inputTreeHits->Add(path_to_hits + "/hits.root");
  std::cout << "Total entries in hits tree: " << inputTreeHits->GetEntries() << std::endl;
  std::cout << "Drawing event_id: " << evToDraw << std::endl;
  int nDrawn;
  if (drawSecondaries)
  {
    std::cout << "Including secondaries in the drawing." << std::endl;
    nDrawn = inputTreeHits->Draw("ty:tx:tz", Form("event_id==%d", evToDraw), "same");
  }
  else
  {
    std::cout << "Drawing only primary hits." << std::endl;
    nDrawn = inputTreeHits->Draw("ty:tx:tz", Form("event_id==%d && barcode_generation==0", evToDraw), "same");
  }
  std::cout << "Number of hits drawn: " << nDrawn << std::endl;

  TGraph *graph_hits = (TGraph *)gPad->GetPrimitive("Graph");
  if (graph_hits != 0x0)
  {
    graph_hits->SetTitle(";x, mm;y, mm");
    graph_hits->SetName(Form("2D_hits"));
    graph_hits->SetMarkerColor(kRed);
    graph_hits->SetMarkerStyle(24);
    graph_hits->SetMarkerSize(1.1);
    leg->AddEntry(graph_hits, "Hits", "p");
  }
  else // 3D
  {
    TPolyMarker3D *graph_hits = (TPolyMarker3D *)gPad->GetPrimitive("TPolyMarker3D");
    drawPM(24, 1.1, kRed, false);
    leg->AddEntry(graph_hits, "Hits", "p");
  }

  // ### states
  if (0)
  {
    inputTreeStates = new TChain("trackstates");
    inputTreeStates->Add(path_to_hits + "/trackstates_ambi.root");

    // loop over "tracks" (for loopers it means "segments")
    int colorsStates[] = {
      kMagenta, kBlue + 1, kOrange + 7, kGray + 2, kCyan, kAzure + 1, kYellow + 1,
      kMagenta, kBlue + 1, kOrange + 7, kGray + 2, kCyan, kAzure + 1, kYellow + 1,
      kMagenta, kBlue + 1, kOrange + 7, kGray + 2, kCyan, kAzure + 1, kYellow + 1,
    };
    for (int k = 0; k < 20; k++)
    {
      inputTreeStates->Draw("t_y:t_x:t_z", Form("layer_id%%2==0 && event_nr==%d && track_nr==%d", evToDraw, k), "same");
      if (k == 0)
      {
          TPolyMarker3D *gr = (TPolyMarker3D *)gPad->GetPrimitive("TPolyMarker3D");
          leg->AddEntry(gr, "track states", "pl");
      }

      drawPM(20, 0.6, colorsStates[k]);
    }
  }

  canv_hits->cd(1);

  gPad->SetGrid();
  leg->Draw();

  // ========== PAD 2: r-z View ==========
  canv_hits->cd(2);
  
  // Draw r vs z where r = sqrt(tx^2 + ty^2)
  int nDrawn_rz, nDrawn_rz_other_event;
  unsigned evOther = 99;
  if (drawSecondaries) {
    nDrawn_rz = inputTreeHits->Draw("sqrt(tx*tx + ty*ty):tz",
                                    Form("event_id==%d", evToDraw), "P");
  } else {
    nDrawn_rz = inputTreeHits->Draw("sqrt(tx*tx + ty*ty):tz",
                                    Form("event_id==%d && barcode_generation==0", evToDraw), "P");
  }
  std::cout << "Number of hits drawn in r-z view: " << nDrawn_rz << std::endl;
  
  TGraph *graph_rz = (TGraph *)gPad->GetPrimitive("Graph");
  if (graph_rz != 0x0)
  {
    graph_rz->SetTitle("r-z view;z (mm);r (mm)");
    graph_rz->SetName("graph_rz");
    graph_rz->SetMarkerColor(kRed);
    graph_rz->SetMarkerStyle(20);
    graph_rz->SetMarkerSize(1.0);
    graph_rz->GetXaxis()->SetTitleOffset(1.2);
    graph_rz->GetYaxis()->SetTitleOffset(1.2);
  }

  // if (drawSecondaries) {
  //   nDrawn_rz_other_event =
  //     inputTreeHits->Draw("sqrt(tx*tx + ty*ty):tz",
  //                         Form("event_id==%d", evOther), "P SAME");
  // } else {
  //   nDrawn_rz_other_event =
  //     inputTreeHits->Draw("sqrt(tx*tx + ty*ty):tz",
  //                         Form("event_id==%d && barcode_generation==0", evOther), "P SAME");
  // }

  // TGraph *graph_rz_ev2 = (TGraph *)gPad->GetPrimitive("Graph");
  // if (graph_rz_ev2 != 0x0) {
  //   graph_rz_ev2->SetName("graph_rz_ev2");
  //   graph_rz_ev2->SetMarkerColor(kBlue);
  //   graph_rz_ev2->SetMarkerStyle(20);
  //   graph_rz_ev2->SetMarkerSize(1.0);
  // }
  
  gPad->SetGrid();
  // gPad->SetLogy(1);
  // gPad->SetLogx(1);

  canv_hits->SaveAs(Form("figures/hits/Event_%d.pdf", evToDraw));

  return;
}
