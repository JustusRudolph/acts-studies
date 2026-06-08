#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>

// helper functions
using Pos = std::pair<float, float>;
using Positions = std::vector<Pos>;
using IDs = std::pair<unsigned, unsigned>;  // bwd/fwd, disc idx
using HitPosMap = std::unordered_map< /* event -> particle idx -> ID -> list of (x,y) in layer */
  unsigned, std::unordered_map<unsigned, std::unordered_map<unsigned, Positions> > >;

unsigned volume_and_layer_id_to_disc_idx(
  unsigned volume_id, unsigned layer_id) {
  if (volume_id == 8) { // OT backward
    if (layer_id <= 4) return 5;  // most "left"most: largest |z|
    else if (layer_id <= 8) return 4;
    else if (layer_id <= 12) return 3;
  } else if (volume_id == 20) { // ML backward
    if (layer_id <= 4) return 2;  // most "left"most: largest |z|
    else if (layer_id <= 8) return 1;
    else if (layer_id <= 12) return 0;
  } else if (volume_id == 22) { // ML forward
    if (layer_id <= 4) return 6;  // most "left"most: smallest |z|
    else if (layer_id <= 8) return 7;
    else if (layer_id <= 12) return 8;
  } else if (volume_id == 25) { // OT forward
    if (layer_id <= 4) return 9;  // most "left"most: smallest |z|
    else if (layer_id <= 8) return 10;
    else if (layer_id <= 12) return 11;
  } else {
    // barrel hit -- ignore
    // return maximum unsigned to indicate invalid
    return std::numeric_limits<unsigned>::max(); // invalid
  }
  return std::numeric_limits<unsigned>::max(); // not a disc hit
}

IDs disc_idx_to_ID(unsigned disc_idx) {
  if (disc_idx <= 5) return {0, disc_idx};  // backward
  else if (disc_idx <= 11) return {1, disc_idx - 6};  // forward
  else {
    std::cerr << "Invalid disc_idx: " << disc_idx << std::endl;
    return {std::numeric_limits<unsigned>::max(),
            std::numeric_limits<unsigned>::max()}; // invalid
  }
}

bool hit_measurement_match(Pos hit_pos, Pos meas_pos,
                          double tol=0.1) {  // need some tolerance in matching: 100um now
  double dx = hit_pos.first - meas_pos.first;
  double dy = hit_pos.second - meas_pos.second;
  return dx*dx < tol*tol && dy*dy < tol*tol;

}

void disc_coverage(
  const std::vector<TString> pathBases = {"geo_staves_pi_1GeV_eta14-20"},
  const bool onSTBC=false,
  const unsigned debugLevel=0)  // 0 nothing, 1 some, 2 many, 3 all debug prints
{
  const TString base             = onSTBC ? "/data/alice/jrudolph/" : "/home/justus/projects/";
  const TString actso2_output_base = base + "alice/ACTSO2/output/";
  const TString geo_studies_base = base + "alice/acts-studies/geo_studies/";

  // 12 disc layers: 6 backward (negative z), 6 forward (positive z)
  // Within each side, layers are ordered by increasing |z|
  constexpr int nDiscs = 6;
  const std::array<TString, nDiscs> bwdLabels = {"bwd_disc_0", "bwd_disc_1", "bwd_disc_2",
                                                  "bwd_disc_3", "bwd_disc_4", "bwd_disc_5"};
  const std::array<TString, nDiscs> fwdLabels = {"fwd_disc_0", "fwd_disc_1", "fwd_disc_2",
                                                  "fwd_disc_3", "fwd_disc_4", "fwd_disc_5"};

  std::array<float, nDiscs> xyMax = {400.0, 400.0, 400.0, 750.0, 750.0, 750.0};  // mm
  const int nXYBins = 500;

  // --- Histograms ---
  // xy hit positions that reach measurements
  std::array<TH2F*, nDiscs> hHitsXY_bwd, hHitsXY_fwd, hMeasXY_bwd, hMeasXY_fwd;
  // xy hit positions that do NOT reach measurements
  std::array<TH2F*, nDiscs> hMissXY_bwd, hMissXY_fwd;
  // efficiency vs phi
  std::array<TEfficiency*, nDiscs> hEffPhi_bwd, hEffPhi_fwd;

  for (int i = 0; i < nDiscs; i++) {
    hHitsXY_bwd[i] = new TH2F("hHitsXY_" + bwdLabels[i],
                               bwdLabels[i] + " full hitmap;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);
    hHitsXY_fwd[i] = new TH2F("hHitsXY_" + fwdLabels[i],
                               fwdLabels[i] + " full hitmap;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);
    hMeasXY_bwd[i] = new TH2F("hMeasXY_" + bwdLabels[i],
                               bwdLabels[i] + " measurement positions;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);
    hMeasXY_fwd[i] = new TH2F("hMeasXY_" + fwdLabels[i],
                               fwdLabels[i] + " measurement positions;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);
    hMissXY_bwd[i] = new TH2F("hMissXY_" + bwdLabels[i],
                               bwdLabels[i] + " hits not becoming measurements;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);
    hMissXY_fwd[i] = new TH2F("hMissXY_" + fwdLabels[i],
                               fwdLabels[i] + " hits not becoming measurements;x (mm);y (mm)",
                               nXYBins, -xyMax[i], xyMax[i], nXYBins, -xyMax[i], xyMax[i]);

    hEffPhi_bwd[i] = new TEfficiency("hEffPhi_" + bwdLabels[i],
                                     bwdLabels[i] + " efficiency vs #phi;#phi (rad);Efficiency",
                                     100, -TMath::Pi(), TMath::Pi());
    hEffPhi_fwd[i] = new TEfficiency("hEffPhi_" + fwdLabels[i],
                                     fwdLabels[i] + " efficiency vs #phi;#phi (rad);Efficiency",
                                     100, -TMath::Pi(), TMath::Pi());

    hHitsXY_bwd[i]->SetStats(false);
    hHitsXY_fwd[i]->SetStats(false);
    hMeasXY_bwd[i]->SetStats(false);
    hMeasXY_fwd[i]->SetStats(false);
    hMissXY_bwd[i]->SetStats(false);
    hMissXY_fwd[i]->SetStats(false);
  }


  // branch variables
  // hit branches
  unsigned hit_eventId, hit_particle_idx, hit_layer_id, hit_volume_id;
  float hit_x, hit_y;
  // measurement branches
  int meas_eventId, meas_layer_id, meas_volume_id;
  /// can be several particles per measurement (clustered hits), with pgun it's rarely not
  /// exactly one, and with primary single pgun sims, it can only be one.
  std::vector<unsigned>* meas_particle_idx = nullptr;
  float meas_true_x, meas_true_y;

  // maps to fill and clear after each file
  HitPosMap hitPosMap;
  HitPosMap measPosMap;


  // --- Event loop over all input files ---
  for (const TString& pathBase : pathBases) {
    std::cout << "Processing " << pathBase << "..." << std::endl;

    const TString hitsFile         = actso2_output_base + pathBase + "/hits.root";
    const TString measurementsFile = actso2_output_base + pathBase + "/measurements.root";

    TFile* fHits = TFile::Open(hitsFile);
    if (!fHits || fHits->IsZombie()) { std::cerr << "Cannot open " << hitsFile << " — skipping" << std::endl; continue; }
    TFile* fMeas = TFile::Open(measurementsFile);
    if (!fMeas || fMeas->IsZombie()) { std::cerr << "Cannot open " << measurementsFile << " — skipping" << std::endl; fHits->Close(); continue; }

    TTree* tHits = (TTree*)fHits->Get("hits");
    if (!tHits) { std::cerr << "hits tree not found in " << hitsFile << " — skipping" << std::endl; fHits->Close(); fMeas->Close(); continue; }
    TTree* tMeas = (TTree*)fMeas->Get("measurements");
    if (!tMeas) { std::cerr << "measurements tree not found in " << measurementsFile << " — skipping" << std::endl; fHits->Close(); fMeas->Close(); continue; }

    // hit branches
    tHits->SetBranchAddress("event_id", &hit_eventId);
    tHits->SetBranchAddress("barcode_particle", &hit_particle_idx);
    tHits->SetBranchAddress("layer_id", &hit_layer_id);
    tHits->SetBranchAddress("volume_id", &hit_volume_id);
    tHits->SetBranchAddress("tx", &hit_x);
    tHits->SetBranchAddress("ty", &hit_y);
    // measurement branches
    tMeas->SetBranchAddress("event_nr", &meas_eventId);
    tMeas->SetBranchAddress("particles_particle", &meas_particle_idx);
    tMeas->SetBranchAddress("layer_id", &meas_layer_id);
    tMeas->SetBranchAddress("volume_id", &meas_volume_id);
    tMeas->SetBranchAddress("true_x", &meas_true_x);
    tMeas->SetBranchAddress("true_y", &meas_true_y);

    // Loop over hits and fill maps
    Long64_t nHits = tHits->GetEntries();
    unsigned nHitsFilled = 0;
    if (debugLevel > 0) std::cout << "DEBUG (1): starting hits loop, nHits="
                         << nHits << std::flush << std::endl;
    for (Long64_t i = 0; i < nHits; i++) {
      if (i % 10000 == 0  && debugLevel > 2)
        std::cout << "DEBUG (3): hits entry " << i << "/"
                  << nHits << std::flush << std::endl;
      else if (i % 100000 == 0 && debugLevel > 1)
        std::cout << "DEBUG (2): hits entry " << i << "/"
                  << nHits << std::flush << std::endl;
      tHits->GetEntry(i);
      unsigned disc_idx = volume_and_layer_id_to_disc_idx(hit_volume_id, hit_layer_id);
      if (disc_idx == std::numeric_limits<unsigned>::max()) continue; // not a disc hit, ignore

      IDs discIDs = disc_idx_to_ID(disc_idx);
      auto [disc_side, layer_in_side] = discIDs;
      if (layer_in_side >= (unsigned)nDiscs) {
        if (debugLevel > 0)
          std::cerr << "DEBUG ERROR: layer_in_side=" << layer_in_side
                    << " out of range at hits entry " << i
                    << " vol=" << hit_volume_id << " layer=" << hit_layer_id
                    << " disc_idx=" << disc_idx << std::flush << std::endl;
        continue;
      }

      // check if there already exist hits before adding another
      auto itDiscIdx = hitPosMap[hit_eventId][hit_particle_idx].find(disc_idx);
      if (itDiscIdx == hitPosMap[hit_eventId][hit_particle_idx].end()) {
        // no hits for this disc yet: add new entry
        hitPosMap[hit_eventId][hit_particle_idx][disc_idx] =
          {std::make_pair(hit_x, hit_y)};
      } else {
        // already hits for this disc: add to vector
        hitPosMap[hit_eventId][hit_particle_idx][disc_idx].
          push_back(std::make_pair(hit_x, hit_y));
      }

      if (disc_side == 0) { // backward
        hHitsXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
      } else if (disc_side == 1) { // forward
        hHitsXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
      }
      nHitsFilled++;
    }
    if (debugLevel > 0)
      std::cout << "DEBUG (0): hits loop done, hitPosMap has "
                << hitPosMap.size() << " events with " << nHitsFilled
                << " hits filled." << std::flush << std::endl;

    // Loop over measurements and fill maps
    Long64_t nMeas = tMeas->GetEntries();
    unsigned nMeasFilled = 0;
    if (debugLevel > 0) std::cout << "DEBUG (0): starting measurements loop, nMeas="
                                  << nMeas << std::flush << std::endl;
    for (Long64_t i = 0; i < nMeas; i++) {
      if (debugLevel > 2)
        std::cout << "DEBUG (3): meas entry " << i << "/"
                  << nMeas << std::flush << std::endl;
      else if (i % 100000 == 0 && debugLevel > 1)
        std::cout << "DEBUG (2): meas entry " << i << "/"
                  << nMeas << std::flush << std::endl;
      tMeas->GetEntry(i);
      unsigned disc_idx = volume_and_layer_id_to_disc_idx(meas_volume_id, meas_layer_id);
      if (meas_particle_idx->size() == 0) continue; // no associated particle, ignore (shouldn't happen, just safety check)
      unsigned part_idx = meas_particle_idx->front();

      if (disc_idx == std::numeric_limits<unsigned>::max()) continue; // not a disc measurement, ignore

      IDs discIDs = disc_idx_to_ID(disc_idx);
      auto [disc_side, layer_in_side] = discIDs;
      if (layer_in_side >= (unsigned)nDiscs) {
        if (debugLevel > 0)
          std::cerr << "DEBUG ERROR: layer_in_side=" << layer_in_side
                    << " out of range at meas entry " << i
                    << " vol=" << meas_volume_id << " layer=" << meas_layer_id
                    << " disc_idx=" << disc_idx << std::flush << std::endl;
        continue;
      }

      // check if there already exist measurements before adding another
      auto itDiscIdx = measPosMap[meas_eventId][part_idx].find(disc_idx);
      if (itDiscIdx == measPosMap[meas_eventId][part_idx].end()) {
        measPosMap[meas_eventId][part_idx][disc_idx] =
          {std::make_pair(meas_true_x, meas_true_y)};
      } else {
        measPosMap[meas_eventId][part_idx][disc_idx].
          push_back(std::make_pair(meas_true_x, meas_true_y));
      }

      if (debugLevel > 2)
        std::cout << "DEBUG (3): meas entry " << i << " eventId=" << meas_eventId
                  << " part_idx=" << part_idx << " disc_idx=" << disc_idx
                  << " meas_true_x=" << meas_true_x << " meas_true_y=" << meas_true_y
                  << std::flush << std::endl;
      if (disc_side == 0) { // backward
        hMeasXY_bwd[layer_in_side]->Fill(meas_true_x, meas_true_y);
      } else if (disc_side == 1) { // forward
        hMeasXY_fwd[layer_in_side]->Fill(meas_true_x, meas_true_y);
      }
      nMeasFilled++;
    }  // loop over measurements

    if (debugLevel > 0)
      std::cout << "DEBUG (0): meas loop done, measPosMap has "
                << measPosMap.size() << " events with " << nMeasFilled
                << " measurements filled." << std::flush << std::endl;

    // loop over events and particles in hitPosMap, check if they have measurement, and fill histograms
    if (debugLevel > 0) std::cout << "DEBUG (0): starting match loop over "
                                  << hitPosMap.size() << " events" << std::flush << std::endl;
    unsigned matchLoopEventCount = 0;
    for (const auto& [eventId, particlesInEvent] : hitPosMap) {
      if (matchLoopEventCount % 1000 == 0 && debugLevel > 2)
        std::cout << "DEBUG (3): match loop event " << matchLoopEventCount
                  << " (eventId=" << eventId << ")" << std::flush << std::endl;
      else if (matchLoopEventCount % 10000 == 0 && debugLevel > 1)
        std::cout << "DEBUG (2): match loop event " << matchLoopEventCount
                  << " (eventId=" << eventId << ")" << std::flush << std::endl;
      matchLoopEventCount++;
      auto itMeasEvent = measPosMap.find(eventId);

      if (itMeasEvent == measPosMap.end()) {
        // ------------------------------- NO MEASUREMENTS IN THIS EVENT -------------------------------
        // no measurement for this event: all hits are misses (this should NEVER happen, but just in case)
        for (const auto& [particleIdx, particleHitLayers] : particlesInEvent) {
          for (const auto& [disc_idx, hit_positions] : particleHitLayers) {
            IDs discIDs = disc_idx_to_ID(disc_idx);
            auto [disc_side, layer_in_side] = discIDs;
            for (const auto& hit_pos : hit_positions) {
              const auto& [hit_x, hit_y] = hit_pos;
              float phi = std::atan2(hit_y, hit_x);
              if (disc_side == 0) { // backward
                hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_bwd[layer_in_side]->Fill(0, phi);
              } else if (disc_side == 1) { // forward
                hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_fwd[layer_in_side]->Fill(0, phi);
              }
            }  // loop over hits in layer 
          }  // loop over layers of particle
        }  // loop over particles in event
        continue;  // move to next event in hitPosMap
      }  // event not found in measurements

      for (const auto& [particleIdx, particleHitLayers] : particlesInEvent) {
        // check if this particle has measurement
        auto itMeasParticle = measPosMap[eventId].find(particleIdx);
        if (itMeasParticle == measPosMap[eventId].end()) {
          // -------------------------------- NO MEASUREMENTS FOR THIS PARTICLE -------------------------------
          // no measurements for this particle: all hits are misses (should almost never happen)
          for (const auto& [disc_idx, hit_positions] : particleHitLayers) {
            IDs discIDs = disc_idx_to_ID(disc_idx);
            auto [disc_side, layer_in_side] = discIDs;
            for (const auto& hit_pos : hit_positions) {
              const auto& [hit_x, hit_y] = hit_pos;
              float phi = std::atan2(hit_y, hit_x);
              if (disc_side == 0) { // backward
                hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_bwd[layer_in_side]->Fill(0, phi);
              } else if (disc_side == 1) { // forward
                hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_fwd[layer_in_side]->Fill(0, phi);
              }
            }  // loop over hits in layer
          }  // loop over layers of particle
          continue;  // move to next particle
        }  // particle not found in measurements
        // -------------------------------- THIS PARTICLE HAS MEASUREMENTS -------------------------------
        auto& measLayers = itMeasParticle->second;
        for (const auto& [disc_idx, hit_positions] : particleHitLayers) {
          auto itMeasDiscIdx = measLayers.find(disc_idx);
          IDs discIDs = disc_idx_to_ID(disc_idx);
          auto [disc_side, layer_in_side] = discIDs;
          if (layer_in_side >= (unsigned)nDiscs) {
            if (debugLevel > 0)
              std::cerr << "DEBUG ERROR: layer_in_side=" << layer_in_side
                        << " out of range in match loop, eventId=" << eventId
                        << " particleIdx=" << particleIdx << " disc_idx=" << disc_idx
                        << std::flush << std::endl;
            continue;
          }
          if (itMeasDiscIdx == measLayers.end()) {
            // -------------------------------- THIS LAYER HAS NO MEASUREMENT -------------------------------
            // no measurement for this layer: all hits are misses (somewhat rare)
            for (const auto& hit_pos : hit_positions) {
              const auto& [hit_x, hit_y] = hit_pos;
              float phi = std::atan2(hit_y, hit_x);
              if (disc_side == 0) { // backward
                hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_bwd[layer_in_side]->Fill(0, phi);
              } else if (disc_side == 1) { // forward
                hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                hEffPhi_fwd[layer_in_side]->Fill(0, phi);
              }
            }  // loop over hits in layer
          }  // no measurement in layer (if statement)
          else {
            // --------------- THIS LAYER HAS AT LEAST ONE MEASUREMENT ----------------------
            // go through each hit to see if it matches any measurement
            auto& meas_positions = itMeasDiscIdx->second;  // can be several per layer
            for (const auto& hit_pos : hit_positions) {
              const auto& [hit_x, hit_y] = hit_pos;
              float phi = std::atan2(hit_y, hit_x);
              bool has_match = false;
              for (auto itMeasPos = meas_positions.begin(); itMeasPos != meas_positions.end(); ++itMeasPos) {
                if (hit_measurement_match(hit_pos, *itMeasPos)) {
                  has_match = true;
                  // remove measurement after match since we do this with single particle event scans
                  meas_positions.erase(itMeasPos);
                  break;
                }
              }
              if (has_match) {
                // hit becomes measurement
                if (disc_side == 0) { // backward
                  hEffPhi_bwd[layer_in_side]->Fill(1, phi);
                } else if (disc_side == 1) { // forward
                  hEffPhi_fwd[layer_in_side]->Fill(1, phi);
                }
              } else {
                // hit does not become measurement: fill miss hist
                if (disc_side == 0) { // backward
                  hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                  hEffPhi_bwd[layer_in_side]->Fill(0, phi);
                } else if (disc_side == 1) { // forward
                  hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                  hEffPhi_fwd[layer_in_side]->Fill(0, phi);
                }
              }
            }  // loop over hits in layer
          }  // layer has at least one measurement (else statement)
        }  // loop over layers of particle
      }  // loop over particles in event
    }  // loop over events
    if (debugLevel > 0) std::cout << "DEBUG (0): match loop done" << std::flush << std::endl;

    fHits->Close();
    fMeas->Close();

    // clear maps for next file
    hitPosMap.clear();
    measPosMap.clear();
    if (debugLevel > 0) std::cout << "DEBUG (0): maps cleared, moving to next file"
                                  << std::flush << std::endl;
  }

  if (debugLevel > 0) std::cout << "DEBUG (0): all files processed, starting drawing"
                                << std::flush << std::endl;

  // --- Draw Canvases ---
  // all on different canvases, sorted by path: efficiencies into figures/geo/efficiency_layer
  // hit maps into figures/geo/occupancy_layer
  TString occupancyOutDir = geo_studies_base + "figures/geo/occupancy_layer/";
  TString efficiencyOutDir = geo_studies_base + "figures/geo/efficiency_layer/";
  for (int i = 0; i < nDiscs; i++) {
    // backward
    TCanvas* c_bwd = new TCanvas("c_bwd_" + bwdLabels[i], bwdLabels[i] + " hit coverage", 1200, 400);
    c_bwd->Divide(3,1);
    c_bwd->cd(1);
    hHitsXY_bwd[i]->Draw("colz");
    c_bwd->cd(2);
    hMeasXY_bwd[i]->Draw("colz");
    c_bwd->cd(3);
    hMissXY_bwd[i]->Draw("colz");
    c_bwd->SaveAs(occupancyOutDir + "xy_hitmap_" + bwdLabels[i] + ".pdf");

    TCanvas* c_eff_bwd = new TCanvas("c_eff_bwd_" + bwdLabels[i], bwdLabels[i] + " efficiency vs phi", 900, 700);
    c_eff_bwd->SetLeftMargin(0.12);
    c_eff_bwd->SetRightMargin(0.15);
    c_eff_bwd->SetTopMargin(0.08);
    c_eff_bwd->SetBottomMargin(0.12);
    hEffPhi_bwd[i]->Draw("AP");
    gPad->Update();
    hEffPhi_bwd[i]->GetPaintedGraph()->GetXaxis()->SetTitle("#phi (rad)");
    hEffPhi_bwd[i]->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
    hEffPhi_bwd[i]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.95, 1.02);
    gPad->Update();
    c_eff_bwd->SaveAs(efficiencyOutDir + "eff_phi_" + bwdLabels[i] + ".pdf");

    // forward
    TCanvas* c_fwd = new TCanvas("c_fwd_" + fwdLabels[i], fwdLabels[i] + " hit coverage", 1200, 400);
    c_fwd->Divide(3,1);
    c_fwd->cd(1);
    hHitsXY_fwd[i]->Draw("colz");
    c_fwd->cd(2);
    hMeasXY_fwd[i]->Draw("colz");
    c_fwd->cd(3);
    hMissXY_fwd[i]->Draw("colz");
    c_fwd->SaveAs(occupancyOutDir + "xy_hitmap_" + fwdLabels[i] + ".pdf");

    TCanvas* c_eff_fwd = new TCanvas("c_eff_fwd_" + fwdLabels[i], fwdLabels[i] + " efficiency vs phi", 900, 700);
    c_eff_fwd->SetLeftMargin(0.12);
    c_eff_fwd->SetRightMargin(0.15);
    c_eff_fwd->SetTopMargin(0.08);
    c_eff_fwd->SetBottomMargin(0.12);
    hEffPhi_fwd[i]->Draw("AP");
    gPad->Update();
    hEffPhi_fwd[i]->GetPaintedGraph()->GetXaxis()->SetTitle("#phi (rad)");
    hEffPhi_fwd[i]->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
    hEffPhi_fwd[i]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.95, 1.02);
    gPad->Update();
    c_eff_fwd->SaveAs(efficiencyOutDir + "eff_phi_" + fwdLabels[i] + ".pdf");
  }
}
