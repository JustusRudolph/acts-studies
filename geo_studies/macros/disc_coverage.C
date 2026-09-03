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
using Pos = std::array<float, 3>;  // x, y, z
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
                           bool include_z = true,
                           double tol=0.1) {  // need some tolerance in matching: 100um now
  double dx = hit_pos[0] - meas_pos[0];
  double dy = hit_pos[1] - meas_pos[1];
  double dz = hit_pos[2] - meas_pos[2];
  bool xy_match = dx*dx < tol*tol && dy*dy < tol*tol;
  if (include_z) {
    return xy_match && dz*dz < tol*tol;
  } else {
    return xy_match;
  }
}

void disc_coverage(
  const std::vector<TString> pathBases = {"geo_staves_pi_1GeV_eta14-20"},
  const bool onSTBC=false,
  const float collision_rate=2400.,  // in kHz (change for PbPb)
  const unsigned nEvents=1000000,
  const bool runWithMeasurements=false,
  const float toleranceML_mm = 3.4,
  const float toleranceOT_mm = 3.4,  // by how much we can go outside nominal outer radius
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
  const std::array<TString, nDiscs> bwdLabelsInPlot =
    {"Backward disc 0 (ML)", "Backward disc 1 (ML)", "Backward disc 2 (ML)",
     "Backward disc 3 (OT)", "Backward disc 4 (OT)", "Backward disc 5 (OT)"};
  const std::array<TString, nDiscs> fwdLabels = {"fwd_disc_0", "fwd_disc_1", "fwd_disc_2",
                                                  "fwd_disc_3", "fwd_disc_4", "fwd_disc_5"};
  const std::array<TString, nDiscs> fwdLabelsInPlot =
    {"Forward disc 0 (ML)", "Forward disc 1 (ML)", "Forward disc 2 (ML)",
     "Forward disc 3 (OT)", "Forward disc 4 (OT)", "Forward disc 5 (OT)"};

  std::array<float, nDiscs> xyMax = {400.0, 400.0, 400.0, 750.0, 750.0, 750.0};  // mm
  std::array<float, nDiscs> rMinNominal = {100.0, 100.0, 100.0, 200.0, 200.0, 200.0};  // mm
  std::array<float, nDiscs> rMaxNominal = {350.0, 350.0, 350.0, 680.0, 680.0, 680.0};  // mm
  std::array<float, nDiscs> rMaxWithTolerance;
  for (int i = 0; i < nDiscs; i++) {
    if (i < 3) { // ML
      rMaxWithTolerance[i] = rMaxNominal[i] + toleranceML_mm;
    } else { // OT
      rMaxWithTolerance[i] = rMaxNominal[i] + toleranceOT_mm;
    }
  }
  // VERY IMPORTANT TO KEEP BIN SIZE TO 1MM FOR COMPARISON WITH STAVE LAYOUTS!
  float binSize_mm = 1.0;  // mm
  std::array<unsigned, nDiscs> nXYBins;
  for (int i = 0; i < nDiscs; i++) {
    nXYBins[i] = static_cast<unsigned>(2 * xyMax[i] / binSize_mm);
  }
  if (debugLevel >= 1) {
    std::cout << "Using bin size of " << binSize_mm << " mm, resulting in nXYBins: ";
    for (int i = 0; i < nDiscs; i++) {
      std::cout << nXYBins[i] << " ";
    }
    std::cout << std::endl;
  }

  // --- Histograms ---
  // xy hit positions that reach measurements
  std::array<TH2F*, nDiscs> hHitsXY_bwd, hHitsXY_fwd, hMeasXY_bwd, hMeasXY_fwd;
  std::array<TH2F*, nDiscs> hHitsXY_primary_bwd, hHitsXY_primary_fwd;
  // xy hit positions that do NOT reach measurements
  std::array<TH2F*, nDiscs> hMissXY_bwd, hMissXY_fwd;
  // nHits wrt x and r
  std::array<TH1D*, nDiscs> hnHits_wrtX_bwd, hnHits_wrtX_fwd,
                            hnHits_wrtX_primary_bwd, hnHits_wrtX_primary_fwd,
                            hnHits_wrtX_Scaled_bwd, hnHits_wrtX_Scaled_fwd,
                            hnHits_wrtX_Scaled_primary_bwd, hnHits_wrtX_Scaled_primary_fwd,
                            hnHits_wrtR_bwd, hnHits_wrtR_fwd,
                            hnHits_wrtR_primary_bwd, hnHits_wrtR_primary_fwd;
  std::array<TEfficiency*, nDiscs> effDisc_bwd_wrt_r, effDisc_fwd_wrt_r;

  TEfficiency* effDisc_bwd_nominal = new TEfficiency(
    "effDisc_bwd", "Hit efficiency by backward layer;Backward disc layer;Efficiency",
    nDiscs, -0.5, nDiscs - 0.5);
  TEfficiency* effDisc_fwd_nominal = new TEfficiency(
    "effDisc_fwd", "Hit efficiency by forward layer;Forward disc layer;Efficiency",
    nDiscs, -0.5, nDiscs - 0.5);
  TEfficiency* effDisc_bwd_outer_ring = new TEfficiency(
    "effDisc_bwd_outer_ring", "Hit efficiency by backward layer (outer ring);Backward disc layer;Efficiency",
    nDiscs, -0.5, nDiscs - 0.5);
  TEfficiency* effDisc_fwd_outer_ring = new TEfficiency(
    "effDisc_fwd_outer_ring", "Hit efficiency by forward layer (outer ring);Forward disc layer;Efficiency",
    nDiscs, -0.5, nDiscs - 0.5);

  for (int i = 0; i < nDiscs; i++) {
    hHitsXY_bwd[i] = new TH2F("hHitsXY_" + bwdLabels[i],
                               bwdLabels[i] + " full hitmap;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hHitsXY_fwd[i] = new TH2F("hHitsXY_" + fwdLabels[i],
                               fwdLabels[i] + " full hitmap;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hMeasXY_bwd[i] = new TH2F("hMeasXY_" + bwdLabels[i],
                               bwdLabels[i] + " measurement positions;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hMeasXY_fwd[i] = new TH2F("hMeasXY_" + fwdLabels[i],
                               fwdLabels[i] + " measurement positions;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hMissXY_bwd[i] = new TH2F("hMissXY_" + bwdLabels[i],
                               bwdLabels[i] + " hits not becoming measurements;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hMissXY_fwd[i] = new TH2F("hMissXY_" + fwdLabels[i],
                               fwdLabels[i] + " hits not becoming measurements;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hHitsXY_primary_bwd[i] = new TH2F("hHitsXY_primary_" + bwdLabels[i],
                               bwdLabels[i] + " hitmap from primaries only;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);
    hHitsXY_primary_fwd[i] = new TH2F("hHitsXY_primary_" + fwdLabels[i],
                               fwdLabels[i] + " hitmap from primaries only;x (mm);y (mm)",
                               nXYBins[i], -xyMax[i], xyMax[i], nXYBins[i], -xyMax[i], xyMax[i]);

    hnHits_wrtX_bwd[i] = new TH1D("hnHits_wrtX_" + bwdLabels[i],
                                  bwdLabels[i] + " number of hits vs x;x (mm);n hits",
                                  nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_primary_bwd[i] = new TH1D("hnHits_wrtX_primary_" + bwdLabels[i],
                                  bwdLabels[i] + " hits from primary particles vs x;x (mm);n hits",
                                  nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_fwd[i] = new TH1D("hnHits_wrtX_" + fwdLabels[i],
                                  fwdLabels[i] + " number of hits vs x;x (mm);n hits",
                                  nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_primary_fwd[i] = new TH1D("hnHits_wrtX_primary_" + fwdLabels[i],
                                  fwdLabels[i] + " hits from primary particles vs x;x (mm);n hits",
                                  nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_Scaled_bwd[i] =
      new TH1D("hnHits_wrtX_Scaled_" + bwdLabels[i],
               bwdLabels[i] + " number of hits vs x;x (mm);n hits",
               nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_Scaled_primary_bwd[i] =
      new TH1D("hnHits_wrtX_Scaled_primary_" + bwdLabels[i],
               bwdLabels[i] + " hits from primary particles vs x;x (mm);n hits",
               nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_Scaled_fwd[i] =
      new TH1D("hnHits_wrtX_Scaled_" + fwdLabels[i],
               fwdLabels[i] + " number of hits vs x;x (mm);n hits",
               nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtX_Scaled_primary_fwd[i] =
      new TH1D("hnHits_wrtX_Scaled_primary_" + fwdLabels[i],
               fwdLabels[i] + " hits from primary particles vs x;x (mm);n hits",
               nXYBins[i], -xyMax[i], xyMax[i]);
    hnHits_wrtR_bwd[i] = new TH1D("hnHits_wrtR_" + bwdLabels[i],
                                  bwdLabels[i] + " number of hits vs r;r (mm);n hits",
                                  nXYBins[i], 0, xyMax[i]);  // not exactly 1mm bins
    hnHits_wrtR_primary_bwd[i] = new TH1D("hnHits_wrtR_primary_" + bwdLabels[i],
                                  bwdLabels[i] + " hits from primary particles vs r;r (mm);n hits",
                                  nXYBins[i], 0, xyMax[i]);  // not exactly 1mm bins
    hnHits_wrtR_fwd[i] = new TH1D("hnHits_wrtR_" + fwdLabels[i],
                                  fwdLabels[i] + " number of hits vs r;r (mm);n hits",
                                  nXYBins[i], 0, xyMax[i]);  // not exactly 1mm bins
    hnHits_wrtR_primary_fwd[i] = new TH1D("hnHits_wrtR_primary_" + fwdLabels[i],
                                  fwdLabels[i] + " hits from primary particles vs r;r (mm);n hits",
                                  nXYBins[i], 0, xyMax[i]);  // not exactly 1mm bins

    effDisc_bwd_wrt_r[i] = new TEfficiency(
      "effDisc_bwd_" + bwdLabels[i],
      bwdLabelsInPlot[i] + " hit efficiency;x (mm);Efficiency",
      50, 0, xyMax[i]);
    effDisc_fwd_wrt_r[i] = new TEfficiency(
      "effDisc_fwd_" + fwdLabels[i],
      fwdLabelsInPlot[i] + " hit efficiency;x (mm);Efficiency",
      50, 0, xyMax[i]);

    hHitsXY_bwd[i]->SetStats(false);
    hHitsXY_fwd[i]->SetStats(false);
    hMeasXY_bwd[i]->SetStats(false);
    hMeasXY_fwd[i]->SetStats(false);
    hMissXY_bwd[i]->SetStats(false);
    hMissXY_fwd[i]->SetStats(false);
    hHitsXY_primary_bwd[i]->SetStats(false);
    hHitsXY_primary_fwd[i]->SetStats(false);

    hnHits_wrtX_bwd[i]->SetStats(false);
    hnHits_wrtX_fwd[i]->SetStats(false);
    hnHits_wrtX_primary_bwd[i]->SetStats(false);
    hnHits_wrtX_primary_fwd[i]->SetStats(false);
    hnHits_wrtX_Scaled_bwd[i]->SetStats(false);
    hnHits_wrtX_Scaled_fwd[i]->SetStats(false);
    hnHits_wrtX_Scaled_primary_bwd[i]->SetStats(false);
    hnHits_wrtX_Scaled_primary_fwd[i]->SetStats(false);
    hnHits_wrtR_bwd[i]->SetStats(false);
    hnHits_wrtR_fwd[i]->SetStats(false);
    hnHits_wrtR_primary_fwd[i]->SetStats(false);
    hnHits_wrtR_primary_bwd[i]->SetStats(false);

  }


  // branch variables
  // hit branches
  unsigned hit_eventId, hit_particle_idx, hit_layer_id, hit_volume_id, hit_part_gen;
  float hit_x, hit_y, hit_z;
  // measurement branches
  int meas_eventId, meas_layer_id, meas_volume_id;
  /// can be several particles per measurement (clustered hits), with pgun it's rarely not
  /// exactly one, and with primary single pgun sims, it can only be one.
  std::vector<unsigned>* meas_particle_idx = nullptr;
  float meas_true_x, meas_true_y, meas_true_z;

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
    tHits->SetBranchAddress("barcode_generation", &hit_part_gen);
    tHits->SetBranchAddress("layer_id", &hit_layer_id);
    tHits->SetBranchAddress("volume_id", &hit_volume_id);
    tHits->SetBranchAddress("tx", &hit_x);
    tHits->SetBranchAddress("ty", &hit_y);
    tHits->SetBranchAddress("tz", &hit_z);
    // measurement branches
    tMeas->SetBranchAddress("event_nr", &meas_eventId);
    tMeas->SetBranchAddress("particles_particle", &meas_particle_idx);
    tMeas->SetBranchAddress("layer_id", &meas_layer_id);
    tMeas->SetBranchAddress("volume_id", &meas_volume_id);
    tMeas->SetBranchAddress("true_x", &meas_true_x);
    tMeas->SetBranchAddress("true_y", &meas_true_y);
    tMeas->SetBranchAddress("true_z", &meas_true_z);

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
          {std::array<float, 3>{hit_x, hit_y, hit_z}};
      } else {
        // already hits for this disc: add to vector
        hitPosMap[hit_eventId][hit_particle_idx][disc_idx].
          push_back(std::array<float, 3>{hit_x, hit_y, hit_z});
      }

      if (disc_side == 0) { // backward
        hHitsXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
        hnHits_wrtX_bwd[layer_in_side]->Fill(hit_x);
        hnHits_wrtR_bwd[layer_in_side]->Fill(std::sqrt(hit_x*hit_x + hit_y*hit_y));
        if (hit_part_gen == 0) { // primary particle
          hHitsXY_primary_bwd[layer_in_side]->Fill(hit_x, hit_y);
          hnHits_wrtX_primary_bwd[layer_in_side]->Fill(hit_x);
          hnHits_wrtR_primary_bwd[layer_in_side]->Fill(std::sqrt(hit_x*hit_x + hit_y*hit_y));
        }
      } else if (disc_side == 1) { // forward
        hHitsXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
        hnHits_wrtX_fwd[layer_in_side]->Fill(hit_x);
        hnHits_wrtR_fwd[layer_in_side]->Fill(std::sqrt(hit_x*hit_x + hit_y*hit_y));
        if (hit_part_gen == 0) { // primary particle
          hHitsXY_primary_fwd[layer_in_side]->Fill(hit_x, hit_y);
          hnHits_wrtX_primary_fwd[layer_in_side]->Fill(hit_x);
          hnHits_wrtR_primary_fwd[layer_in_side]->Fill(std::sqrt(hit_x*hit_x + hit_y*hit_y));
        }
      }
      nHitsFilled++;
    }
    if (debugLevel > 0)
      std::cout << "DEBUG (0): hits loop done, hitPosMap has "
                << hitPosMap.size() << " events with " << nHitsFilled
                << " hits filled." << std::flush << std::endl;

    if (runWithMeasurements) {  // comparison to measurements only if asked for
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
            {std::array<float, 3>{meas_true_x, meas_true_y, meas_true_z}};
        } else {
          measPosMap[meas_eventId][part_idx][disc_idx].
            push_back(std::array<float, 3>{meas_true_x, meas_true_y, meas_true_z});
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
                const auto& [hit_x, hit_y, hit_z] = hit_pos;
                float phi = std::atan2(hit_y, hit_x);
                if (disc_side == 0) { // backward
                  hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                } else if (disc_side == 1) { // forward
                  hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
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
                const auto& [hit_x, hit_y, hit_z] = hit_pos;
                float phi = std::atan2(hit_y, hit_x);
                if (disc_side == 0) { // backward
                  hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                } else if (disc_side == 1) { // forward
                  hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
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
                const auto& [hit_x, hit_y, hit_z] = hit_pos;
                float phi = std::atan2(hit_y, hit_x);
                if (disc_side == 0) { // backward
                  hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                } else if (disc_side == 1) { // forward
                  hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                }
              }  // loop over hits in layer
            }  // no measurement in layer (if statement)
            else {
              // --------------- THIS LAYER HAS AT LEAST ONE MEASUREMENT ----------------------
              // go through each hit to see if it matches any measurement
              auto& meas_positions = itMeasDiscIdx->second;  // can be several per layer
              for (const auto& hit_pos : hit_positions) {
                const auto& [hit_x, hit_y, hit_z] = hit_pos;
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
                  } else if (disc_side == 1) { // forward
                  }
                } else {
                  // hit does not become measurement: fill miss hist
                  if (disc_side == 0) { // backward
                    hMissXY_bwd[layer_in_side]->Fill(hit_x, hit_y);
                  } else if (disc_side == 1) { // forward
                    hMissXY_fwd[layer_in_side]->Fill(hit_x, hit_y);
                  }
                }
              }  // loop over hits in layer
            }  // layer has at least one measurement (else statement)
          }  // loop over layers of particle
        }  // loop over particles in event
      }  // loop over events
      if (debugLevel > 0) std::cout << "DEBUG (0): match loop done" << std::flush << std::endl;
    }  // if runWithMeasurements

    fHits->Close();
    fMeas->Close();

    // clear maps for next file
    hitPosMap.clear();
    measPosMap.clear();
    if (debugLevel > 0) std::cout << "DEBUG (0): maps cleared, moving to next file"
                                  << std::flush << std::endl;
  }
  // now time to scale with r and stave length
  // an excellent proxy for stave length is the number of non-zero bins in the y
  // projection at a given x bin. This is the same nominally for all ML and OT respectively
  std::vector<double> staveLengthsX_ML(nXYBins[0]), staveLengthsX_OT(nXYBins[3]);
  // use most filled histos for stave length "calculation"
  for (unsigned i_xy_bin_ml = 0; i_xy_bin_ml < nXYBins[0]; i_xy_bin_ml++) {
    // get y projection of ML histos at this x bin and count non-empty bins
    staveLengthsX_ML[i_xy_bin_ml] = 0;
    for (int i_y_bin = 1; i_y_bin <= nXYBins[0]; i_y_bin++) {
      if (hnHits_wrtX_bwd[0]->GetBinContent(i_xy_bin_ml, i_y_bin) > 0)
        staveLengthsX_ML[i_xy_bin_ml]++;
    }
  } // loop over x bins for ML
  for (unsigned i_xy_bin_ot = 0; i_xy_bin_ot < nXYBins[3]; i_xy_bin_ot++) {
    // get y projection of OT histos at this x bin and count non-empty bins
    staveLengthsX_OT[i_xy_bin_ot] = 0;
    for (int i_y_bin = 1; i_y_bin <= nXYBins[3]; i_y_bin++) {
      if (hnHits_wrtX_fwd[3]->GetBinContent(i_xy_bin_ot, i_y_bin) > 0)
        staveLengthsX_OT[i_xy_bin_ot]++;
    }
  } // loop over x bins for OT

  // perform fits on r dependence after all hits are written to later scale x dep
  std::array<TF1*, nDiscs> fitFunc_bwd, fitFunc_fwd;
  std::array<TF1*, nDiscs> fitFunc_primary_bwd, fitFunc_primary_fwd;
  for (int i = 0; i < nDiscs; i++) {
    // make fit in reduced range to avoid edge effects
    fitFunc_bwd[i] = new TF1("fitFunc_bwd_" + bwdLabels[i], "[0]/pow(x, [1])",
                              rMinNominal[i] + 50, rMaxNominal[i] - 50);
    // start with 1/r dependence
    fitFunc_bwd[i]->SetParameters(hnHits_wrtR_bwd[i]->GetMean(), 1);
    // r for restricting range to given r values above
    hnHits_wrtR_bwd[i]->Fit(fitFunc_bwd[i], "QR");

    fitFunc_primary_bwd[i] = new TF1("fitFunc_primary_bwd_" + bwdLabels[i], "[0]/pow(x, [1])",
                              rMinNominal[i] + 50, rMaxNominal[i] - 50);
    fitFunc_primary_bwd[i]->SetParameters(hnHits_wrtR_primary_bwd[i]->GetMean(), 1);
    hnHits_wrtR_primary_bwd[i]->Fit(fitFunc_primary_bwd[i], "QR");

    fitFunc_fwd[i] = new TF1("fitFunc_fwd_" + fwdLabels[i], "[0]/pow(x, [1])",
                              rMinNominal[i] + 50, rMaxNominal[i] - 50);
    // start with 1/r dependence
    fitFunc_fwd[i]->SetParameters(hnHits_wrtR_fwd[i]->GetMean(), 1);
    hnHits_wrtR_fwd[i]->Fit(fitFunc_fwd[i], "QR");

    fitFunc_primary_fwd[i] = new TF1("fitFunc_primary_fwd_" + fwdLabels[i], "[0]/pow(x, [1])",
                              rMinNominal[i] + 50, rMaxNominal[i] - 50);
    fitFunc_primary_fwd[i]->SetParameters(hnHits_wrtR_primary_fwd[i]->GetMean(), 1);
    hnHits_wrtR_primary_fwd[i]->Fit(fitFunc_primary_fwd[i], "QR");
  }

  // loop over hits again for scaled x dependence
  for (const TString& pathBase : pathBases) {
    const TString hitsFile         = actso2_output_base + pathBase + "/hits.root";
    TFile* fHits = TFile::Open(hitsFile);
    if (!fHits || fHits->IsZombie()) {
      std::cerr << "Cannot open " << hitsFile << " for second loop — skipping" << std::endl;
      continue;
    }
    TTree* tHits = (TTree*)fHits->Get("hits");
    if (!tHits) {
      std::cerr << "hits tree not found in " << hitsFile
                << " for second loop — skipping" << std::endl; fHits->Close();
      continue;
    }
    // hit branches
    tHits->SetBranchAddress("layer_id", &hit_layer_id);
    tHits->SetBranchAddress("volume_id", &hit_volume_id);
    tHits->SetBranchAddress("tx", &hit_x);
    tHits->SetBranchAddress("ty", &hit_y);

    unsigned nHits = tHits->GetEntries();
    if (debugLevel > 0) std::cout << "DEBUG (1): hits loop again after fitting: nHits="
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
      double r = std::sqrt(hit_x*hit_x + hit_y*hit_y);
      double scale_factor_stave_length = 1.0;
      if (layer_in_side < 3) { // ML
        scale_factor_stave_length = staveLengthsX_ML[hnHits_wrtX_bwd[0]->FindBin(hit_x)];
      } else { // OT
        scale_factor_stave_length = staveLengthsX_OT[hnHits_wrtX_fwd[3]->FindBin(hit_x)];
      }
      // "convert" from mm to m, effectively just give neater numbers in plot
      scale_factor_stave_length /= 1000;
      if (disc_side == 0) { // backward
        double scaleFactor = fitFunc_bwd[layer_in_side]->Eval(r);
        hnHits_wrtX_Scaled_bwd[layer_in_side]->Fill(
          hit_x, 1.0/(scaleFactor * scale_factor_stave_length));
        if (hit_part_gen == 0) { // primary particle
          double scaleFactor = fitFunc_primary_bwd[layer_in_side]->Eval(r);
          hnHits_wrtX_Scaled_primary_bwd[layer_in_side]->Fill(
            hit_x, 1.0/(scaleFactor * scale_factor_stave_length));
        }
      } else if (disc_side == 1) { // forward
        double scaleFactor = fitFunc_fwd[layer_in_side]->Eval(r);
        hnHits_wrtX_Scaled_fwd[layer_in_side]->Fill(
          hit_x, 1.0/(scaleFactor * scale_factor_stave_length));
        if (hit_part_gen == 0) { // primary particle
          double scaleFactor = fitFunc_primary_fwd[layer_in_side]->Eval(r);
          hnHits_wrtX_Scaled_primary_fwd[layer_in_side]->Fill(
            hit_x, 1.0/(scaleFactor * scale_factor_stave_length));
        }
      }
    }  // loop over hits again for scaled x dependence
    fHits->Close();
  }  // loop over files again for scaled x dependence

  // now fill efficiency by layer histograms after we have all the hits processed
  for (unsigned i_disc = 0; i_disc < nDiscs; i_disc++) {
    for (unsigned i_xbin = 0; i_xbin < nXYBins[i_disc]; i_xbin++) {
      for (unsigned i_ybin = 0; i_ybin < nXYBins[i_disc]; i_ybin++) {
        unsigned nHits_bwd = hHitsXY_bwd[i_disc]->GetBinContent(i_xbin, i_ybin);
        unsigned nHits_fwd = hHitsXY_fwd[i_disc]->GetBinContent(i_xbin, i_ybin);
        // bin positions are identical for backward and forward since they have the same binning
        float x_mm = hHitsXY_bwd[i_disc]->GetXaxis()->GetBinCenter(i_xbin);
        float y_mm = hHitsXY_bwd[i_disc]->GetYaxis()->GetBinCenter(i_ybin);
        float r_mm = std::sqrt(x_mm * x_mm + y_mm * y_mm);
        effDisc_bwd_wrt_r[i_disc]->Fill(nHits_bwd > 0, r_mm);
        effDisc_fwd_wrt_r[i_disc]->Fill(nHits_fwd > 0, r_mm);
        if (r_mm < rMinNominal[i_disc]) {
          continue;
        } else if (r_mm < rMaxNominal[i_disc]) {
          // larger than inner, smaller than outer: fill nominal eff histos
          effDisc_bwd_nominal->Fill(nHits_bwd > 0, i_disc);
          effDisc_fwd_nominal->Fill(nHits_fwd > 0, i_disc);
        } else if (r_mm < rMaxWithTolerance[i_disc]) {
          // larger than outer, smaller than extended/with tolerance: fill extended eff histos
          effDisc_bwd_outer_ring->Fill(nHits_bwd > 0, i_disc);
          effDisc_fwd_outer_ring->Fill(nHits_fwd > 0, i_disc);
        }  // region checking
      }  // loop over y bins
    }  // loop over x bins
  }  // loop over discs

  if (debugLevel > 0) std::cout << "DEBUG (0): all files processed, starting drawing"
                                << std::flush << std::endl;

  // --- Draw Canvases ---
  // all on different canvases, sorted by path: efficiencies into figures/geo/efficiency_layer
  // hit maps into figures/geo/occupancy_layer
  // for everything, only make colz log when there are at least 100 files
  TString occupancyOutDir = geo_studies_base + "figures/geo/occupancy_layer/";
  TString efficiencyOutDir = geo_studies_base + "figures/geo/efficiency_layer/";
  for (int i = 0; i < nDiscs; i++) {
    // backward
    TCanvas* c_bwd;
    if (runWithMeasurements) {
      c_bwd = new TCanvas("c_bwd_" + bwdLabels[i], bwdLabels[i] + " hit coverage", 1400, 400);
      c_bwd->Divide(3,1);
      c_bwd->cd(1);
      hHitsXY_bwd[i]->SetMinimum(1); // set minimum for log scale to 1
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hHitsXY_bwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hHitsXY_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": Hit positions;x (mm);y (mm)");
      c_bwd->cd(2);
      hMeasXY_bwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hMeasXY_bwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hMeasXY_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": Measurement positions;x (mm);y (mm)");
      c_bwd->cd(3);
      hMissXY_bwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hMissXY_bwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hMissXY_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": Hits not becoming measurements;x (mm);y (mm)");
    } else {
      c_bwd = new TCanvas("c_bwd_" + bwdLabels[i], bwdLabels[i] + " hit coverage", 600, 500);
      hHitsXY_bwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hHitsXY_bwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hHitsXY_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": Hit positions;x (mm);y (mm)");
    }
    c_bwd->SaveAs(occupancyOutDir + "xy_hitmap_" + bwdLabels[i] + ".pdf");

    TCanvas* c_bwd_primary = new TCanvas(
      "c_bwd_primary_" + bwdLabels[i], bwdLabels[i] + " hit coverage from primaries", 600, 500);
    hHitsXY_primary_bwd[i]->SetMinimum(1);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
    hHitsXY_primary_bwd[i]->Draw("colz");
    if (pathBases.size() >= 100) gPad->SetLogz();
    hHitsXY_primary_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": Hit positions from primaries;x (mm);y (mm)");
    c_bwd_primary->SaveAs(occupancyOutDir + "xy_hitmap_primary_" + bwdLabels[i] + ".pdf");

    TCanvas* c_hits1D_bwd = new TCanvas("c_hits1D_bwd_" + bwdLabels[i], bwdLabels[i] + " nHits wrt x and r", 1400, 900);
    // split into three: r dep hits + fit, x dep hits, x dep scaled hits
    c_hits1D_bwd->Divide(3,2);
    c_hits1D_bwd->cd(1);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtR_bwd[i]->Draw("EP");
    hnHits_wrtR_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": N_{hits} vs r;r (mm);N_{hits}");
    fitFunc_bwd[i]->Draw("SAME");
    // legend for fit
    TLegend* legend_bwd = new TLegend(0.6, 0.6, 0.85, 0.85);
    legend_bwd->AddEntry(hnHits_wrtR_bwd[i], "Hits", "lep");
    TString fitLabel = Form("Fit: #frac{%.1f}{r^{%.2f}}",
                            fitFunc_bwd[i]->GetParameter(0),
                            fitFunc_bwd[i]->GetParameter(1));
    legend_bwd->AddEntry(fitFunc_bwd[i], fitLabel, "l");
    legend_bwd->SetBorderSize(0);
    legend_bwd->SetFillStyle(0);
    legend_bwd->Draw();

    c_hits1D_bwd->cd(2);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_bwd[i]->Draw("EP");
    hnHits_wrtX_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": N_{hits} vs x;x (mm);N_{hits}");

    c_hits1D_bwd->cd(3);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_Scaled_bwd[i]->Draw("EP");
    hnHits_wrtX_Scaled_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": N_{hits} vs x (scaled);x (mm);Relative hit density");

    // now the same but for primaries only
    c_hits1D_bwd->cd(4);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtR_primary_bwd[i]->Draw("EP");
    hnHits_wrtR_primary_bwd[i]->SetTitle(bwdLabelsInPlot[i] + ": hits from primaries vs r;r (mm);N_{hits}");
    fitFunc_primary_bwd[i]->Draw("SAME");
    // legend for fit
    TLegend* legend_primary_bwd = new TLegend(0.6, 0.6, 0.85, 0.85);
    legend_primary_bwd->AddEntry(hnHits_wrtR_primary_bwd[i], "Hits", "lep");
    TString fitLabel_primary_bwd = Form("Fit: #frac{%.1f}{r^{%.2f}}",
                                        fitFunc_primary_bwd[i]->GetParameter(0),
                                        fitFunc_primary_bwd[i]->GetParameter(1));
    legend_primary_bwd->AddEntry(fitFunc_primary_bwd[i], fitLabel_primary_bwd, "l");
    legend_primary_bwd->SetBorderSize(0);
    legend_primary_bwd->SetFillStyle(0);
    legend_primary_bwd->Draw();

    c_hits1D_bwd->cd(5);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_primary_bwd[i]->Draw("EP");
    std::cout << "number of entries in primary bwd hits vs x for " << bwdLabels[i] << ": "
              << hnHits_wrtX_primary_bwd[i]->GetEntries() << std::endl;
    hnHits_wrtX_primary_bwd[i]->SetTitle(
      bwdLabelsInPlot[i] + ": hits from primaries vs x;x (mm);N_{hits}");

    c_hits1D_bwd->cd(6);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_Scaled_primary_bwd[i]->Draw("EP");
    hnHits_wrtX_Scaled_primary_bwd[i]->SetTitle(
      bwdLabelsInPlot[i] + ": hits from primaries vs x (scaled);x (mm);Relative hit density");

    c_hits1D_bwd->SaveAs(efficiencyOutDir + "nHits_wrt_X_R_" + bwdLabels[i] + ".pdf");

    TCanvas* c_eff_bwd = new TCanvas("c_eff_bwd_" + bwdLabels[i], bwdLabels[i] + " efficiency vs r", 600, 500);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    effDisc_bwd_wrt_r[i]->Draw("EP");
    effDisc_bwd_wrt_r[i]->SetTitle(bwdLabelsInPlot[i] + ": Efficiency vs r;r (mm);Efficiency");
    c_eff_bwd->SaveAs(efficiencyOutDir + "efficiency_vs_r_" + bwdLabels[i] + ".pdf");

    // forward
    TCanvas* c_fwd;
    if (runWithMeasurements) {
      c_fwd = new TCanvas("c_fwd_" + fwdLabels[i], fwdLabels[i] + " hit coverage", 1400, 400);
      c_fwd->Divide(3,1);
      c_fwd->cd(1);
      hHitsXY_fwd[i]->SetMinimum(1); // set minimum for log scale to 1
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hHitsXY_fwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hHitsXY_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": Hit positions;x (mm);y (mm)");
      c_fwd->cd(2);
      hMeasXY_fwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hMeasXY_fwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hMeasXY_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": Measurement positions;x (mm);y (mm)");
      c_fwd->cd(3);
      hMissXY_fwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hMissXY_fwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hMissXY_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": Hits not becoming measurements;x (mm);y (mm)");
    } else {
      c_fwd = new TCanvas("c_fwd_" + fwdLabels[i], fwdLabels[i] + " hit coverage", 600, 500);
      hHitsXY_fwd[i]->SetMinimum(1);
      gPad->SetLeftMargin(0.15);  // y label cut off otherwise
      gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
      hHitsXY_fwd[i]->Draw("colz");
      if (pathBases.size() >= 100) gPad->SetLogz();
      hHitsXY_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": Hit positions;x (mm);y (mm)");
    }
    c_fwd->SaveAs(occupancyOutDir + "xy_hitmap_" + fwdLabels[i] + ".pdf");

    TCanvas* c_fwd_primary = new TCanvas(
      "c_fwd_primary_" + fwdLabels[i], fwdLabels[i] + " hit coverage for primaries", 600, 500);
    hHitsXY_primary_fwd[i]->SetMinimum(1);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    gPad->SetRightMargin(0.15);  // important for colz — leaves room for the z-axis palette
    hHitsXY_primary_fwd[i]->Draw("colz");
    if (pathBases.size() >= 100) gPad->SetLogz();
    hHitsXY_primary_fwd[i]->SetTitle(
      fwdLabelsInPlot[i] + ": Hit positions from primaries only;x (mm);y (mm)");
    c_fwd_primary->SaveAs(occupancyOutDir + "xy_hitmap_primaries_" + fwdLabels[i] + ".pdf");

    TCanvas* c_hits1D_fwd = new TCanvas(
      "c_hits1D_fwd_" + fwdLabels[i], fwdLabels[i] + " nHits wrt x and r", 1400, 900);
    // split into three: r dep hits + fit, x dep hits, x dep scaled hits
    c_hits1D_fwd->Divide(3,2);
    c_hits1D_fwd->cd(1);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtR_fwd[i]->Draw("EP");
    hnHits_wrtR_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": N_{hits} vs r;r (mm);N_{hits}");
    fitFunc_fwd[i]->Draw("SAME");
    // legend for fit
    TLegend* legend_fwd = new TLegend(0.6, 0.6, 0.85, 0.85);
    legend_fwd->AddEntry(hnHits_wrtR_fwd[i], "Hits", "lep");
    TString fitLabel_fwd = Form("Fit: #frac{%.1f}{r^{%.2f}}",
                                fitFunc_fwd[i]->GetParameter(0),
                                fitFunc_fwd[i]->GetParameter(1));
    legend_fwd->AddEntry(fitFunc_fwd[i], fitLabel_fwd, "l");
    legend_fwd->SetBorderSize(0);
    legend_fwd->SetFillStyle(0);
    legend_fwd->Draw();

    c_hits1D_fwd->cd(2);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_fwd[i]->Draw("EP");
    hnHits_wrtX_fwd[i]->SetTitle(fwdLabelsInPlot[i] + ": N_{hits} vs x;x (mm);N_{hits}");

    c_hits1D_fwd->cd(3);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_Scaled_fwd[i]->Draw("EP");
    hnHits_wrtX_Scaled_fwd[i]->SetTitle(
      fwdLabelsInPlot[i] + ": N_{hits} vs x (scaled);x (mm);Relative hit density");
    
    c_hits1D_fwd->cd(4);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtR_primary_fwd[i]->Draw("EP");
    hnHits_wrtR_primary_fwd[i]->SetTitle(
      fwdLabelsInPlot[i] + ": hits from primaries vs r;r (mm);N_{hits}");
    fitFunc_primary_fwd[i]->Draw("SAME");
    // legend for fit
    TLegend* legend_primary_fwd = new TLegend(0.6, 0.6, 0.85, 0.85);
    legend_primary_fwd->AddEntry(hnHits_wrtR_primary_fwd[i], "Hits", "lep");
    TString fitLabel_primary_fwd = Form("Fit: #frac{%.1f}{r^{%.2f}}",
                                        fitFunc_primary_fwd[i]->GetParameter(0),
                                        fitFunc_primary_fwd[i]->GetParameter(1));
    legend_primary_fwd->AddEntry(fitFunc_primary_fwd[i], fitLabel_primary_fwd, "l");
    legend_primary_fwd->SetBorderSize(0);
    legend_primary_fwd->SetFillStyle(0);
    legend_primary_fwd->Draw();

    c_hits1D_fwd->cd(5);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_primary_fwd[i]->Draw("EP");
    hnHits_wrtX_primary_fwd[i]->SetTitle(
      fwdLabelsInPlot[i] + ": hits from primaries vs x;x (mm);N_{hits}");

    c_hits1D_fwd->cd(6);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    hnHits_wrtX_Scaled_primary_fwd[i]->Draw("EP");
    hnHits_wrtX_Scaled_primary_fwd[i]->SetTitle(
      fwdLabelsInPlot[i] + ": hits from primaries vs x (scaled);x (mm);Relative hit density");
    
    c_hits1D_fwd->SaveAs(efficiencyOutDir + "nHits_wrt_X_R_" + fwdLabels[i] + ".pdf");

    TCanvas* c_eff_fwd = new TCanvas("c_eff_fwd_" + fwdLabels[i], fwdLabels[i] + " efficiency vs r", 600, 500);
    gPad->SetLeftMargin(0.15);  // y label cut off otherwise
    effDisc_fwd_wrt_r[i]->Draw("EP");
    effDisc_fwd_wrt_r[i]->SetTitle(fwdLabelsInPlot[i] + ": Efficiency vs r;r (mm);Efficiency");
    c_eff_fwd->SaveAs(efficiencyOutDir + "efficiency_vs_r_" + fwdLabels[i] + ".pdf");
  }
  // now draw the layer hit efficiency histograms
  TCanvas* c_eff = new TCanvas("c_eff", "Hit efficiency by layer", 1400, 600);
  c_eff->Divide(2,1);
  c_eff->cd(1);
  effDisc_bwd_nominal->Draw("EP");
  effDisc_bwd_nominal->SetLineColor(kBlue);
  effDisc_bwd_nominal->SetTitle("Hit efficiency by layer for nominal acceptance;\
                                 Layer;Efficiency");
  effDisc_fwd_nominal->Draw("EP SAME");
  effDisc_fwd_nominal->SetLineColor(kRed);
  // legend for both
  TLegend* legend_hiteff = new TLegend(0.6, 0.6, 0.85, 0.85);
  legend_hiteff->AddEntry(effDisc_bwd_nominal, "Backward", "lep");
  legend_hiteff->AddEntry(effDisc_fwd_nominal, "Forward", "lep");
  legend_hiteff->SetBorderSize(0);
  legend_hiteff->SetFillStyle(0);
  legend_hiteff->Draw();
  c_eff->cd(2);
  effDisc_bwd_outer_ring->Draw("EP");
  effDisc_bwd_outer_ring->SetLineColor(kBlue);
  effDisc_bwd_outer_ring->SetTitle("Hit efficiency by layer for tolerance region;\
                                    Layer;Efficiency");
  effDisc_fwd_outer_ring->Draw("EP SAME");
  effDisc_fwd_outer_ring->SetLineColor(kRed);
  legend_hiteff->Draw();
  c_eff->SaveAs(efficiencyOutDir + "hit_efficiency_by_layer.pdf");
}
