#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../../Utils/BasicFunctions.h"

// middle layers used for the doublet study, besides the inner TOF; ordered by radius
constexpr unsigned kMiddleVolume = 40;
const std::vector<unsigned> kMiddleLayers = {9, 11, 13, 15, 19};  // r = 70, 90, 120, 190, 300 mm
constexpr unsigned kMinMiddleLayers = 3;

// outer tracker barrel, layers 3, 5, 7 at r = 442, 617, 792 mm
constexpr unsigned kOuterTrackerVolume = 50;

constexpr int kPhotonPdg = 22;
constexpr unsigned kPhotonConversionProcess = 2;  // ActsFatras::GenerationProcess::ePhotonConversion

namespace {

// event, vertex_primary, vertex_secondary, particle, generation, sub_particle
using ParticleKey = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;

// a simulated hit
struct SimHit {
  int index;  // position of the hit along the particle trajectory
  unsigned volume;
  unsigned layer;
  float x;
  float y;
  float z;
  float pt;  // true pT at the hit [GeV]
  double r() const { return std::hypot(x, y); }
  double phi() const { return std::atan2(y, x); }
};

// a selected conversion e+/e- and its hits, ordered along the trajectory
struct ConversionParticle {
  int pdg;
  float pt;
  std::vector<SimHit> hits;
};

bool isSameLayer(const SimHit& a, const SimHit& b) {
  return a.volume == b.volume && a.layer == b.layer;
}

// cot(theta) = dz/dr of the doublet of two hits (symmetric in the two)
double doubletCotTheta(const SimHit& a, const SimHit& b) {
  return (b.z - a.z) / (b.r() - a.r());
}

// cot(theta) of the top doublet minus that of the bottom doublet, triplet ordered bottom, middle, top
double tripletDeltaCotTheta(const std::array<SimHit, 3>& triplet) {
  return doubletCotTheta(triplet[1], triplet[2]) - doubletCotTheta(triplet[0], triplet[1]);
}

// radial position of a middle layer in kMiddleLayers, kMiddleLayers.size() if not a middle layer
int middleLayerRank(const SimHit& hit) {
  if (hit.volume != kMiddleVolume) return kMiddleLayers.size();
  return std::find(kMiddleLayers.begin(), kMiddleLayers.end(), hit.layer) - kMiddleLayers.begin();
}

// middle layer hits of the first outward leg: along the trajectory, one hit per
// layer, each in a layer further out than the previous, until the particle turns
// back inward. Conversions at larger r often go out through the outer layers, curl
// and come back in through the inner ones; mixing both legs gives no valid seed.
std::vector<SimHit> outgoingMiddleLayerHits(const std::vector<SimHit>& hits) {
  std::vector<SimHit> outgoing_hits;
  for (const SimHit& hit : hits) {
    const int rank = middleLayerRank(hit);
    if (rank == int(kMiddleLayers.size())) continue;  // not a middle layer
    if (outgoing_hits.empty() || rank > middleLayerRank(outgoing_hits.back())) {
      outgoing_hits.push_back(hit);
    } else if (hit.layer != outgoing_hits.back().layer) {
      break;  // turned back inward
    }  // else: second hit in the same layer (module overlap), skip
  }
  return outgoing_hits;
}

// whether the particle crosses at least kMinMiddleLayers middle layers on its first outward leg
bool passesMiddleLayers(const std::vector<SimHit>& hits) {
  return outgoingMiddleLayerHits(hits).size() >= kMinMiddleLayers;
}

// phi difference wrapped to [-pi, pi)
double deltaPhi(double phi_to, double phi_from) {
  return std::remainder(phi_to - phi_from, 2 * M_PI);
}

// detector region of a set of hits: all in the middle layers (incl. the inner TOF),
// all in the outer tracker, or anything else (both, or other layers/volumes)
enum Region { kRegionMiddle = 0, kRegionOuter = 1, kRegionOther = 2 };
const char* kRegionLabels[] = {"ML", "OT", "Other"};

Region hitRegion(const SimHit& hit, unsigned tof_volume, unsigned tof_layer) {
  const bool isTof = hit.volume == tof_volume && hit.layer == tof_layer;
  if (isTof || middleLayerRank(hit) < int(kMiddleLayers.size())) return kRegionMiddle;
  if (hit.volume == kOuterTrackerVolume) return kRegionOuter;
  return kRegionOther;
}

template <size_t N>
bool hasOtherLayerHit(const std::array<SimHit, N>& hits, unsigned tof_volume, unsigned tof_layer) {
  return std::any_of(hits.begin(), hits.end(), [&](const SimHit& hit) {
    return hitRegion(hit, tof_volume, tof_layer) == kRegionOther;
  });
}

// mean true pT [GeV] at the hits of a doublet/triplet
template <size_t N>
double meanHitPt(const std::array<SimHit, N>& hits) {
  double sum = 0;
  for (const SimHit& hit : hits) sum += hit.pt;
  return sum / N;
}

template <size_t N>
Region hitsRegion(const std::array<SimHit, N>& hits, unsigned tof_volume, unsigned tof_layer) {
  const Region first_region = hitRegion(hits[0], tof_volume, tof_layer);
  for (const SimHit& hit : hits) {
    if (hitRegion(hit, tof_volume, tof_layer) != first_region) return kRegionOther;
  }
  return first_region;
}

// consecutive hits along the trajectory in pairwise different layers: starting at
// each hit, take the next N; if all in different layers use them and jump past
// them, otherwise move on by one hit. Each group is ordered by radius like a seed.
template <size_t N>
std::vector<std::array<SimHit, N>> consecutiveHitGroups(const std::vector<SimHit>& hits) {
  std::vector<std::array<SimHit, N>> groups;
  size_t i_hit = 0;
  while (i_hit + N <= hits.size()) {
    bool isAllDifferentLayers = true;
    for (size_t i_first = i_hit; i_first < i_hit + N; i_first++) {
      for (size_t i_second = i_first + 1; i_second < i_hit + N; i_second++) {
        if (isSameLayer(hits[i_first], hits[i_second])) isAllDifferentLayers = false;
      }
    }
    if (!isAllDifferentLayers) {
      i_hit++;
      continue;
    }
    std::array<SimHit, N> group;
    std::copy(hits.begin() + i_hit, hits.begin() + i_hit + N, group.begin());
    std::sort(group.begin(), group.end(), [](const SimHit& a, const SimHit& b) { return a.r() < b.r(); });
    groups.push_back(group);
    i_hit += N;
  }
  return groups;
}

// circle through the transverse positions of three hits
struct Circle {
  double x_center;
  double y_center;
  double radius;
  double rCenter() const { return std::hypot(x_center, y_center); }
};

std::optional<Circle> circleThroughHits(const std::array<SimHit, 3>& triplet) {
  const double ax = triplet[0].x, ay = triplet[0].y;
  const double bx = triplet[1].x, by = triplet[1].y;
  const double cx = triplet[2].x, cy = triplet[2].y;
  const double denominator = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  if (std::abs(denominator) < 1e-6) return std::nullopt;  // collinear
  const double a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
  const double x_center = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / denominator;
  const double y_center = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / denominator;
  return Circle{x_center, y_center, std::hypot(ax - x_center, ay - y_center)};
}

// pT [GeV] of a circle of radius [mm] in a field [T]
double circlePt(double radius, double b_field) {
  return 0.299792458e-3 * b_field * radius;
}

// the triplets before and after a turn of a looper, each ordered by radius
// (bottom, middle, top) like a seed, so both use the same layers in the same roles
struct LooperTriplets {
  std::array<SimHit, 3> before_turn;
  std::array<SimHit, 3> after_turn;
};

// first turn along the trajectory: hits n and n+1 in the same layer, as well as
// n-1 and n+2, and n-2 and n+3, with three different layers in each triplet
std::optional<LooperTriplets> findLooperTurn(const std::vector<SimHit>& hits) {
  for (size_t i_hit = 2; i_hit + 3 < hits.size(); i_hit++) {
    const bool isTurn = isSameLayer(hits[i_hit],     hits[i_hit + 1]) &&
                        isSameLayer(hits[i_hit - 1], hits[i_hit + 2]) &&
                        isSameLayer(hits[i_hit - 2], hits[i_hit + 3]);
    const bool isThreeLayers = !isSameLayer(hits[i_hit - 2], hits[i_hit - 1]) &&
                               !isSameLayer(hits[i_hit - 1], hits[i_hit]) &&
                               !isSameLayer(hits[i_hit - 2], hits[i_hit]);
    if (!isTurn || !isThreeLayers) continue;
    LooperTriplets triplets{{hits[i_hit - 2], hits[i_hit - 1], hits[i_hit]},
                            {hits[i_hit + 3], hits[i_hit + 2], hits[i_hit + 1]}};
    for (std::array<SimHit, 3>* triplet : {&triplets.before_turn, &triplets.after_turn}) {
      std::sort(triplet->begin(), triplet->end(), [](const SimHit& a, const SimHit& b) { return a.r() < b.r(); });
    }
    return triplets;
  }
  return std::nullopt;
}

}  // namespace

// pT distribution, production vertex r and production vertex (r, z) of the
// e+/e- from photon conversions in <input_dir>/particles_simulation.root,
// input_dir being relative to the ACTSO2 output base. Additionally the pT of
// those that leave a hit in the inner TOF, identified by its (volume, layer) id
// in <input_dir>/hits.root.
//
// For the conversion e+/e- with an inner TOF hit and hits in at least
// kMinMiddleLayers of the middle layers on their first outward leg, the doublets
// from their true hits: phi difference of the 2nd middle layer hit to the 1st
// and to the 3rd. Seeds may also start at a later hit, so this is only an
// indication of the doublet phi windows needed. Likewise the difference in
// cot(theta) = dz/dr between the top (2nd -> 3rd hit) and bottom (1st -> 2nd hit) doublet.
//
// All consecutive hit triplets (doublets) of all conversion e+/e- with three (two)
// different layers, see consecutiveHitGroups: cot(theta) difference of the
// triplets and phi difference (outer - inner) of the doublets, separately for
// triplets/doublets fully in the middle layers (ML, incl. inner TOF), fully in the
// outer tracker barrel (OT), and all others (Other: ML + OT, or any hit in another
// layer, e.g. the vertex detector or the discs). OT ones only if the true pT at each
// of their hits is at least ot_pt_min. Drawn for all pT and split into bins of the
// mean true pT at their hits, edges group_pt_edges.
//
// Loopers, among all conversion e+/e-: at the first turn where the particle hits
// the same three layers in reverse order (see findLooperTurn), the triplet before
// and after the turn are compared: the sum of their cot(theta) differences (which
// cancel for an ideal helix, both triplets ordered by radius like a seed), and the
// difference in pT estimate and in radial position of the centre of the circle
// through their transverse positions, (after - before) / before.
//
// Selected are generation 1 e+/e- whose parent is a photon and whose creation
// process is a photon conversion. The parent is not stored in the particle
// file; for generation 1 it is the generation 0 particle with the same
// (vertex_primary, vertex_secondary, particle) barcode. The process cut on top
// of the photon parent removes the Compton/photoelectric e- of a photon.
//
// Usage (from low_pt_studies):
//   root -l -b -q 'macros/gen1_electron_vertices.C("pythia_2500ev/seed_0")'
void gen1_electron_vertices(const TString input_dir="pythia_2500ev/seed_0",
                            const bool onSTBC=false,
                            // inner TOF (volume id, layer id), currently at r = 211 mm
                            const std::tuple<unsigned, unsigned> tof_volume_layer={40, 17},
                            const double pt_min=1e-3,   // GeV, log binning
                            const double pt_max=10.,    // GeV
                            const double r_max=1500.,   // mm
                            const double z_max=4000.,   // mm, |z|
                            const unsigned n_bins_pt=100,
                            const unsigned n_bins_r=300,
                            const unsigned n_bins_z=400,
                            const double delta_phi_max=0.5,  // rad
                            const unsigned n_bins_delta_phi=200,
                            const double delta_cot_theta_max=0.5,
                            const unsigned n_bins_delta_cot_theta=200,
                            const double b_field=2.,  // T, for the looper pT estimate
                            const double looper_sum_delta_cot_theta_max=0.5,
                            const double looper_rel_delta_pt_max=1.,
                            const double looper_rel_delta_center_r_max=1.,
                            const unsigned n_bins_looper=200,
                            const double ot_pt_min=0.14,  // GeV, true pT at the hits of OT doublets/triplets
                            // GeV, bins of the mean true pT at the hits of the consecutive doublets/triplets
                            const std::vector<double> group_pt_edges={0., 0.05, 0.1, 0.14, 0.2, 0.3, 0.5, 1e3}) {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);  // keep the histograms alive after closing the input file
  const auto [tof_volume, tof_layer] = tof_volume_layer;

  // logarithmic pT binning
  std::vector<double> pt_edges(n_bins_pt + 1);
  for (unsigned i_bin = 0; i_bin <= n_bins_pt; i_bin++) {
    pt_edges[i_bin] = pt_min * std::pow(pt_max / pt_min, double(i_bin) / n_bins_pt);
  }
  TH1D* pt_electrons_hist = new TH1D("pt_electrons_hist",
                                     "p_{T} of e^{#pm} from #gamma conversions;p_{T} (GeV/c);Entries",
                                     n_bins_pt, pt_edges.data());
  TH1D* pt_positrons_hist = (TH1D*)pt_electrons_hist->Clone("pt_positrons_hist");
  TH1D* pt_tof_electrons_hist = (TH1D*)pt_electrons_hist->Clone("pt_tof_electrons_hist");
  pt_tof_electrons_hist->SetTitle("p_{T} of e^{#pm} from #gamma conversions with an inner TOF hit");
  TH1D* pt_tof_positrons_hist = (TH1D*)pt_tof_electrons_hist->Clone("pt_tof_positrons_hist");
  TH1D* vertex_r_hist = new TH1D("vertex_r_hist",
                                 "Production vertex r of e^{#pm} from #gamma conversions;r (mm);Entries",
                                 n_bins_r, 0, r_max);
  TH2D* vertex_rz_hist = new TH2D("vertex_rz_hist",
                                  "Production vertex of e^{#pm} from #gamma conversions;z (mm);r (mm);Entries",
                                  n_bins_z, -z_max, z_max, n_bins_r, 0, r_max);
  TH1D* delta_phi_bottom_hist = new TH1D("delta_phi_bottom_hist",
                                         "#Delta#phi of the true middle layer hits, conversion e^{#pm} with an inner TOF hit;"
                                         "#Delta#phi (rad);Entries",
                                         n_bins_delta_phi, -delta_phi_max, delta_phi_max);
  TH1D* delta_phi_top_hist = (TH1D*)delta_phi_bottom_hist->Clone("delta_phi_top_hist");
  TH1D* delta_cot_theta_hist = new TH1D("delta_cot_theta_hist",
                                        "Doublet cot#theta difference, conversion e^{#pm} with an inner TOF hit;"
                                        "cot#theta_{top} - cot#theta_{bottom};Entries",
                                        n_bins_delta_cot_theta, -delta_cot_theta_max, delta_cot_theta_max);
  std::array<TH1D*, 3> hit_doublet_delta_phi_hists;
  std::array<TH1D*, 3> hit_triplet_delta_cot_theta_hists;
  for (int i_region = 0; i_region < 3; i_region++) {
    hit_doublet_delta_phi_hists[i_region] = new TH1D(
        Form("hit_doublet_delta_phi_%s_hist", kRegionLabels[i_region]),
        "#Delta#phi of consecutive hit doublets, conversion e^{#pm};#phi_{outer} - #phi_{inner} (rad);Entries",
        n_bins_delta_phi, -delta_phi_max, delta_phi_max);
    hit_triplet_delta_cot_theta_hists[i_region] = new TH1D(
        Form("hit_triplet_delta_cot_theta_%s_hist", kRegionLabels[i_region]),
        "Doublet cot#theta difference of consecutive hit triplets, conversion e^{#pm};"
        "cot#theta_{top} - cot#theta_{bottom};Entries",
        n_bins_delta_cot_theta, -delta_cot_theta_max, delta_cot_theta_max);
  }
  // the same per bin of the mean true pT at the hits
  const size_t n_group_pt_bins = group_pt_edges.size() - 1;
  std::vector<std::array<TH1D*, 3>> pt_bin_doublet_delta_phi_hists(n_group_pt_bins);
  std::vector<std::array<TH1D*, 3>> pt_bin_triplet_delta_cot_theta_hists(n_group_pt_bins);
  for (size_t i_pt_bin = 0; i_pt_bin < n_group_pt_bins; i_pt_bin++) {
    for (int i_region = 0; i_region < 3; i_region++) {
      pt_bin_doublet_delta_phi_hists[i_pt_bin][i_region] = (TH1D*)hit_doublet_delta_phi_hists[i_region]->Clone(
          Form("hit_doublet_delta_phi_%s_pt%zu_hist", kRegionLabels[i_region], i_pt_bin));
      pt_bin_triplet_delta_cot_theta_hists[i_pt_bin][i_region] = (TH1D*)hit_triplet_delta_cot_theta_hists[i_region]->Clone(
          Form("hit_triplet_delta_cot_theta_%s_pt%zu_hist", kRegionLabels[i_region], i_pt_bin));
    }
  }
  // bin of the mean true pT, n_group_pt_bins if outside the edges
  auto groupPtBin = [&](double group_pt) -> size_t {
    if (group_pt < group_pt_edges.front()) return n_group_pt_bins;
    return std::upper_bound(group_pt_edges.begin(), group_pt_edges.end(), group_pt) - group_pt_edges.begin() - 1;
  };
  unsigned long n_other_triplets_outside_ml_ot = 0;
  TH1D* looper_delta_cot_theta_before_hist = new TH1D("looper_delta_cot_theta_before_hist",
                                                      "Looper triplets before and after the turn, conversion e^{#pm};"
                                                      "#Deltacot#theta = cot#theta_{top} - cot#theta_{bottom};Entries",
                                                      n_bins_looper, -looper_sum_delta_cot_theta_max,
                                                      looper_sum_delta_cot_theta_max);
  TH1D* looper_sum_delta_cot_theta_hist = (TH1D*)looper_delta_cot_theta_before_hist->Clone("looper_sum_delta_cot_theta_hist");
  TH1D* looper_delta_pt_hist = new TH1D("looper_delta_pt_hist",
                                        "Looper triplets, conversion e^{#pm};"
                                        "(p_{T, after turn} - p_{T, before turn}) / p_{T, before turn};Entries",
                                        n_bins_looper, -looper_rel_delta_pt_max, looper_rel_delta_pt_max);
  TH1D* looper_delta_center_r_hist = new TH1D("looper_delta_center_r_hist",
                                              "Looper triplets, conversion e^{#pm};"
                                              "(r_{centre, after turn} - r_{centre, before turn}) / r_{centre, before turn};Entries",
                                              n_bins_looper, -looper_rel_delta_center_r_max, looper_rel_delta_center_r_max);

  // --- conversion e+/e- from the particles -----------------------------------
  const TString sim_filename = Utils::outputPathBase(onSTBC) + input_dir + "/particles_simulation.root";
  TFile* sim_file = TFile::Open(sim_filename);
  if (!sim_file || sim_file->IsZombie()) {
    std::cout << "Error: Could not open " << sim_filename << std::endl;
    return;
  }
  TTree* sim_tree = (TTree*)sim_file->Get("particles");
  if (!sim_tree) {
    std::cout << "Error: particles tree not found in " << sim_filename << std::endl;
    sim_file->Close();
    return;
  }

  unsigned event_id;
  std::vector<unsigned>* vertex_primary = nullptr;
  std::vector<unsigned>* vertex_secondary = nullptr;
  std::vector<unsigned>* particle = nullptr;
  std::vector<unsigned>* sub_particle = nullptr;
  std::vector<int>* particle_type = nullptr;
  std::vector<unsigned>* generation = nullptr;
  std::vector<unsigned>* process = nullptr;
  std::vector<float>* pt = nullptr;
  std::vector<float>* vx = nullptr;
  std::vector<float>* vy = nullptr;
  std::vector<float>* vz = nullptr;
  sim_tree->SetBranchStatus("*", false);
  for (const char* name : {"event_id", "vertex_primary", "vertex_secondary", "particle", "sub_particle",
                           "particle_type", "generation", "process", "pt", "vx", "vy", "vz"}) {
    sim_tree->SetBranchStatus(name, true);
  }
  sim_tree->SetBranchAddress("event_id",         &event_id);
  sim_tree->SetBranchAddress("vertex_primary",   &vertex_primary);
  sim_tree->SetBranchAddress("vertex_secondary", &vertex_secondary);
  sim_tree->SetBranchAddress("particle",         &particle);
  sim_tree->SetBranchAddress("sub_particle",     &sub_particle);
  sim_tree->SetBranchAddress("particle_type",    &particle_type);
  sim_tree->SetBranchAddress("generation",       &generation);
  sim_tree->SetBranchAddress("process",          &process);
  sim_tree->SetBranchAddress("pt",               &pt);
  sim_tree->SetBranchAddress("vx",               &vx);
  sim_tree->SetBranchAddress("vy",               &vy);
  sim_tree->SetBranchAddress("vz",               &vz);

  std::map<ParticleKey, ConversionParticle> conversion_particles;
  unsigned long n_parent_missing = 0;
  const Long64_t n_events = sim_tree->GetEntries();
  for (Long64_t i_ev = 0; i_ev < n_events; i_ev++) {
    sim_tree->GetEntry(i_ev);

    // PDG of the generation 0 particles, by (vertex_primary, vertex_secondary, particle)
    std::map<std::tuple<unsigned, unsigned, unsigned>, int> gen0_pdg;
    for (size_t i_mcp = 0; i_mcp < generation->size(); i_mcp++) {
      if (generation->at(i_mcp) != 0) continue;
      gen0_pdg[{vertex_primary->at(i_mcp), vertex_secondary->at(i_mcp), particle->at(i_mcp)}] =
          particle_type->at(i_mcp);
    }

    for (size_t i_mcp = 0; i_mcp < generation->size(); i_mcp++) {
      if (generation->at(i_mcp) != 1) continue;
      const int pdg = particle_type->at(i_mcp);
      if (std::abs(pdg) != 11) continue;
      if (process->at(i_mcp) != kPhotonConversionProcess) continue;
      const auto parent = gen0_pdg.find({vertex_primary->at(i_mcp), vertex_secondary->at(i_mcp), particle->at(i_mcp)});
      if (parent == gen0_pdg.end()) {
        n_parent_missing++;
        continue;
      }
      if (parent->second != kPhotonPdg) continue;
      (pdg == 11 ? pt_electrons_hist : pt_positrons_hist)->Fill(pt->at(i_mcp));
      const double vertex_r = std::hypot(vx->at(i_mcp), vy->at(i_mcp));
      vertex_r_hist->Fill(vertex_r);
      vertex_rz_hist->Fill(vz->at(i_mcp), vertex_r);
      const ParticleKey key{event_id, vertex_primary->at(i_mcp), vertex_secondary->at(i_mcp),
                            particle->at(i_mcp), generation->at(i_mcp), sub_particle->at(i_mcp)};
      conversion_particles[key] = {pdg, pt->at(i_mcp), {}};
    }
  }  // loop over events
  sim_file->Close();
  delete sim_file;
  if (n_parent_missing) {
    std::cout << "Warning: no parent found for " << n_parent_missing << " conversion e+/e-" << std::endl;
  }

  // --- hits of the conversion e+/e- ------------------------------------------
  {
    const TString hits_filename = Utils::outputPathBase(onSTBC) + input_dir + "/hits.root";
    TFile* hits_file = TFile::Open(hits_filename);
    if (!hits_file || hits_file->IsZombie()) {
      std::cout << "Error: Could not open " << hits_filename << std::endl;
      return;
    }
    TTree* hits_tree = (TTree*)hits_file->Get("hits");
    unsigned hit_event_id, hit_vertex_primary, hit_vertex_secondary, hit_particle, hit_generation, hit_sub_particle;
    unsigned volume_id, layer_id;
    int index;
    float tx, ty, tz, tpx, tpy;
    hits_tree->SetBranchStatus("*", false);
    for (const char* name : {"event_id", "barcode_vertex_primary", "barcode_vertex_secondary",
                             "barcode_particle", "barcode_generation", "barcode_sub_particle",
                             "volume_id", "layer_id", "index", "tx", "ty", "tz", "tpx", "tpy"}) {
      hits_tree->SetBranchStatus(name, true);
    }
    hits_tree->SetBranchAddress("event_id",                 &hit_event_id);
    hits_tree->SetBranchAddress("barcode_vertex_primary",   &hit_vertex_primary);
    hits_tree->SetBranchAddress("barcode_vertex_secondary", &hit_vertex_secondary);
    hits_tree->SetBranchAddress("barcode_particle",         &hit_particle);
    hits_tree->SetBranchAddress("barcode_generation",       &hit_generation);
    hits_tree->SetBranchAddress("barcode_sub_particle",     &hit_sub_particle);
    hits_tree->SetBranchAddress("volume_id",                &volume_id);
    hits_tree->SetBranchAddress("layer_id",                 &layer_id);
    hits_tree->SetBranchAddress("index",                    &index);
    hits_tree->SetBranchAddress("tx",                       &tx);
    hits_tree->SetBranchAddress("ty",                       &ty);
    hits_tree->SetBranchAddress("tz",                       &tz);
    hits_tree->SetBranchAddress("tpx",                      &tpx);
    hits_tree->SetBranchAddress("tpy",                      &tpy);
    for (Long64_t i_hit = 0; i_hit < hits_tree->GetEntries(); i_hit++) {
      hits_tree->GetEntry(i_hit);
      const auto conversion_particle = conversion_particles.find(
          {hit_event_id, hit_vertex_primary, hit_vertex_secondary, hit_particle, hit_generation, hit_sub_particle});
      if (conversion_particle == conversion_particles.end()) continue;
      conversion_particle->second.hits.push_back({index, volume_id, layer_id, tx, ty, tz, float(std::hypot(tpx, tpy))});
    }
    hits_file->Close();
    delete hits_file;
  }

  // --- per conversion e+/e-: inner TOF, doublets, loopers ---------------------
  unsigned long n_doublet_particles = 0;
  unsigned long n_loopers = 0;
  for (auto& [key, conversion_particle] : conversion_particles) {
    std::vector<SimHit>& hits = conversion_particle.hits;
    std::sort(hits.begin(), hits.end(), [](const SimHit& a, const SimHit& b) { return a.index < b.index; });

    // all consecutive hit doublets and triplets in different layers
    auto isBelowOtPtMin = [&](const auto& group) {
      return std::any_of(group.begin(), group.end(), [&](const SimHit& hit) { return hit.pt < ot_pt_min; });
    };
    for (const std::array<SimHit, 2>& doublet : consecutiveHitGroups<2>(hits)) {
      const Region region = hitsRegion(doublet, tof_volume, tof_layer);
      if (region == kRegionOuter && isBelowOtPtMin(doublet)) continue;
      const double delta_phi = deltaPhi(doublet[1].phi(), doublet[0].phi());
      hit_doublet_delta_phi_hists[region]->Fill(delta_phi);
      const size_t i_pt_bin = groupPtBin(meanHitPt(doublet));
      if (i_pt_bin < n_group_pt_bins) pt_bin_doublet_delta_phi_hists[i_pt_bin][region]->Fill(delta_phi);
    }
    for (const std::array<SimHit, 3>& triplet : consecutiveHitGroups<3>(hits)) {
      const Region region = hitsRegion(triplet, tof_volume, tof_layer);
      if (region == kRegionOuter && isBelowOtPtMin(triplet)) continue;
      const double delta_cot_theta = tripletDeltaCotTheta(triplet);
      hit_triplet_delta_cot_theta_hists[region]->Fill(delta_cot_theta);
      const size_t i_pt_bin = groupPtBin(meanHitPt(triplet));
      if (i_pt_bin < n_group_pt_bins) pt_bin_triplet_delta_cot_theta_hists[i_pt_bin][region]->Fill(delta_cot_theta);
      if (region == kRegionOther && hasOtherLayerHit(triplet, tof_volume, tof_layer)) n_other_triplets_outside_ml_ot++;
    }

    // loopers
    if (const std::optional<LooperTriplets> triplets = findLooperTurn(hits)) {
      const std::optional<Circle> circle_before = circleThroughHits(triplets->before_turn);
      const std::optional<Circle> circle_after = circleThroughHits(triplets->after_turn);
      if (circle_before && circle_after) {
        const double delta_cot_theta_before = tripletDeltaCotTheta(triplets->before_turn);
        looper_delta_cot_theta_before_hist->Fill(delta_cot_theta_before);
        looper_sum_delta_cot_theta_hist->Fill(delta_cot_theta_before + tripletDeltaCotTheta(triplets->after_turn));
        const double pt_before = circlePt(circle_before->radius, b_field);
        looper_delta_pt_hist->Fill((circlePt(circle_after->radius, b_field) - pt_before) / pt_before);
        looper_delta_center_r_hist->Fill((circle_after->rCenter() - circle_before->rCenter()) / circle_before->rCenter());
        n_loopers++;
      }
    }

    // inner TOF
    const bool hasTofHit = std::any_of(hits.begin(), hits.end(), [&](const SimHit& hit) {
      return hit.volume == tof_volume && hit.layer == tof_layer;
    });
    if (!hasTofHit) continue;
    (conversion_particle.pdg == 11 ? pt_tof_electrons_hist : pt_tof_positrons_hist)->Fill(conversion_particle.pt);

    // doublets from the true hits in the middle layers
    if (!passesMiddleLayers(hits)) continue;
    const std::vector<SimHit> first_hits = outgoingMiddleLayerHits(hits);
    delta_phi_bottom_hist->Fill(deltaPhi(first_hits[1].phi(), first_hits[0].phi()));
    delta_phi_top_hist->Fill(deltaPhi(first_hits[2].phi(), first_hits[1].phi()));
    delta_cot_theta_hist->Fill(doubletCotTheta(first_hits[1], first_hits[2]) -
                               doubletCotTheta(first_hits[0], first_hits[1]));
    n_doublet_particles++;
  }
  std::cout << "Selected " << conversion_particles.size() << " e+/e- from photon conversions in "
            << n_events << " events ("
            << pt_electrons_hist->GetEntries() << " e-, "
            << pt_positrons_hist->GetEntries() << " e+), of which "
            << pt_tof_electrons_hist->GetEntries() << " e- and "
            << pt_tof_positrons_hist->GetEntries() << " e+ hit the inner TOF (volume "
            << tof_volume << ", layer " << tof_layer << "); "
            << n_doublet_particles << " of those have hits in >= " << kMinMiddleLayers
            << " middle layers on the way out" << std::endl;
  std::cout << "Other consecutive triplets: " << n_other_triplets_outside_ml_ot << " of "
            << hit_triplet_delta_cot_theta_hists[kRegionOther]->GetEntries()
            << " have a hit outside ML and OT (vertex detector, discs, ...)" << std::endl;
  std::cout << n_loopers << " conversion e+/e- loop back through the same three layers" << std::endl;
  for (const TH1D* hist : {delta_phi_bottom_hist, delta_phi_top_hist, delta_cot_theta_hist,
                           looper_delta_cot_theta_before_hist, looper_sum_delta_cot_theta_hist,
                           looper_delta_pt_hist, looper_delta_center_r_hist,
                           hit_doublet_delta_phi_hists[0], hit_doublet_delta_phi_hists[1], hit_doublet_delta_phi_hists[2],
                           hit_triplet_delta_cot_theta_hists[0], hit_triplet_delta_cot_theta_hists[1],
                           hit_triplet_delta_cot_theta_hists[2]}) {
    const double n_outside = hist->GetBinContent(0) + hist->GetBinContent(hist->GetNbinsX() + 1);
    if (n_outside > 0) {
      std::cout << hist->GetName() << ": " << n_outside << " of " << hist->GetEntries()
                << " entries outside the plotted range" << std::endl;
    }
  }

  // figures/gamma_conversions/<input_dir with "/" replaced by "_">
  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/low_pt_studies/figures/gamma_conversions/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);

  // --- pT: e+ + e-, e-, e+ ---------------------------------------------------
  TH1D* pt_all_hist = (TH1D*)pt_electrons_hist->Clone("pt_all_hist");
  pt_all_hist->Add(pt_positrons_hist);

  TCanvas* pt_canvas = new TCanvas("c_gen1_electron_pt", "Conversion e+/e- pT", 800, 600);
  pt_all_hist->SetLineColor(kBlack);
  pt_electrons_hist->SetLineColor(kBlue);
  pt_positrons_hist->SetLineColor(kRed);
  pt_all_hist->SetMinimum(0.5);
  pt_all_hist->Draw("HIST");
  pt_electrons_hist->Draw("HIST SAME");
  pt_positrons_hist->Draw("HIST SAME");
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TLegend* legend = new TLegend(0.70, 0.75, 0.88, 0.88);
  legend->AddEntry(pt_all_hist,       "e^{+} + e^{-}", "l");
  legend->AddEntry(pt_electrons_hist, "e^{-}",         "l");
  legend->AddEntry(pt_positrons_hist, "e^{+}",         "l");
  legend->Draw();
  pt_canvas->SaveAs(Form("%s/gen1_electron_pt.pdf", output_dir.Data()));

  // --- pT with an inner TOF hit: e+ + e-, e-, e+ ----------------------------
  TH1D* pt_tof_all_hist = (TH1D*)pt_tof_electrons_hist->Clone("pt_tof_all_hist");
  pt_tof_all_hist->Add(pt_tof_positrons_hist);

  TCanvas* pt_tof_canvas = new TCanvas("c_gen1_electron_pt_tof", "Conversion e+/e- pT, inner TOF hit", 800, 600);
  pt_tof_all_hist->SetLineColor(kBlack);
  pt_tof_electrons_hist->SetLineColor(kBlue);
  pt_tof_positrons_hist->SetLineColor(kRed);
  pt_tof_all_hist->SetMinimum(0.5);
  pt_tof_all_hist->Draw("HIST");
  pt_tof_electrons_hist->Draw("HIST SAME");
  pt_tof_positrons_hist->Draw("HIST SAME");
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  // fraction of all conversion e+/e- (no pT range restriction) with an inner TOF hit
  auto tofFraction = [](const TH1D* tof_hist, const TH1D* all_hist) {
    return all_hist->GetEntries() > 0 ? 100. * tof_hist->GetEntries() / all_hist->GetEntries() : 0.;
  };
  TLegend* tof_legend = new TLegend(0.55, 0.72, 0.88, 0.88);
  tof_legend->SetHeader("fraction with inner TOF hit");
  tof_legend->AddEntry(pt_tof_all_hist,
                       Form("e^{+} + e^{-}: %.3g%%", tofFraction(pt_tof_all_hist, pt_all_hist)), "l");
  tof_legend->AddEntry(pt_tof_electrons_hist,
                       Form("e^{-}: %.3g%%", tofFraction(pt_tof_electrons_hist, pt_electrons_hist)), "l");
  tof_legend->AddEntry(pt_tof_positrons_hist,
                       Form("e^{+}: %.3g%%", tofFraction(pt_tof_positrons_hist, pt_positrons_hist)), "l");
  tof_legend->Draw();
  pt_tof_canvas->SaveAs(Form("%s/gen1_electron_pt_inner_tof.pdf", output_dir.Data()));

  // --- doublet phi differences ----------------------------------------------
  TCanvas* delta_phi_canvas = new TCanvas("c_gen1_electron_doublet_delta_phi", "Conversion e+/e- doublet delta phi", 800, 600);
  delta_phi_bottom_hist->SetLineColor(kBlue);
  delta_phi_top_hist->SetLineColor(kRed);
  delta_phi_bottom_hist->SetMaximum(1.1 * std::max(delta_phi_bottom_hist->GetMaximum(), delta_phi_top_hist->GetMaximum()));
  delta_phi_bottom_hist->Draw("HIST");
  delta_phi_top_hist->Draw("HIST SAME");
  TLegend* delta_phi_legend = new TLegend(0.60, 0.75, 0.88, 0.88);
  delta_phi_legend->SetTextSize(0.04);
  delta_phi_legend->AddEntry(delta_phi_bottom_hist, "#phi_{2} - #phi_{1}", "l");
  delta_phi_legend->AddEntry(delta_phi_top_hist,    "#phi_{3} - #phi_{2}", "l");
  delta_phi_legend->Draw();
  delta_phi_canvas->SaveAs(Form("%s/gen1_electron_doublet_delta_phi.pdf", output_dir.Data()));

  // --- doublet cot(theta) difference ----------------------------------------
  TCanvas* delta_cot_theta_canvas = new TCanvas("c_gen1_electron_doublet_delta_cot_theta",
                                                "Conversion e+/e- doublet delta cot theta", 800, 600);
  delta_cot_theta_hist->SetLineColor(kBlue);
  delta_cot_theta_hist->Draw("HIST");
  delta_cot_theta_canvas->SaveAs(Form("%s/gen1_electron_doublet_delta_cot_theta.pdf", output_dir.Data()));

  // --- all consecutive hit doublets / triplets: ML, OT, Other on the same axes,
  // one pad for all pT and one per bin of the mean true pT at the hits ---------
  auto ptBinLabel = [&](size_t i_pt_bin) -> TString {
    const double low = 1e3 * group_pt_edges[i_pt_bin], high = 1e3 * group_pt_edges[i_pt_bin + 1];
    if (i_pt_bin == 0 && low <= 0) return Form("#LTp_{T}#GT < %.0f MeV/c", high);
    if (i_pt_bin + 1 == n_group_pt_bins) return Form("#LTp_{T}#GT #geq %.0f MeV/c", low);
    return Form("%.0f #leq #LTp_{T}#GT < %.0f MeV/c", low, high);
  };
  auto drawRegionHists = [&](const std::array<TH1D*, 3>& hists, const TString& pad_label) {
    const Color_t colors[] = {kBlue, kRed, kBlack};
    double maximum = 0;
    for (const TH1D* hist : hists) maximum = std::max(maximum, hist->GetMaximum());
    TLegend* region_legend = new TLegend(0.50, 0.70, 0.89, 0.89);
    region_legend->SetBorderSize(0);
    region_legend->SetFillStyle(0);
    for (int i_region = 0; i_region < 3; i_region++) {
      TH1D* hist = hists[i_region];
      hist->SetLineColor(colors[i_region]);
      hist->SetMinimum(0.5);
      hist->SetMaximum(std::max(2 * maximum, 1.));
      hist->Draw(i_region == 0 ? "HIST" : "HIST SAME");
      const TString label = i_region == kRegionOuter ? Form("OT, p_{T} #geq %.0f MeV/c", 1e3 * ot_pt_min)
                                                     : TString(kRegionLabels[i_region]);
      const double n_outside = hist->GetBinContent(0) + hist->GetBinContent(hist->GetNbinsX() + 1);
      const TString count = n_outside > 0 ? Form("%.0f, %.0f%% outside", hist->GetEntries(),
                                                 100. * n_outside / hist->GetEntries())
                                          : Form("%.0f", hist->GetEntries());
      region_legend->AddEntry(hist, Form("%s (%s)", label.Data(), count.Data()), "l");
    }
    region_legend->Draw();
    TPaveText* pad_title = new TPaveText(0.12, 0.80, 0.48, 0.89, "NDC");
    pad_title->SetBorderSize(0);
    pad_title->SetFillStyle(0);
    pad_title->SetTextAlign(12);
    pad_title->AddText(pad_label);
    pad_title->Draw();
    gPad->SetLogy(1);
  };
  auto drawPtBinnedRegionHists = [&](const std::array<TH1D*, 3>& all_pt_hists,
                                     const std::vector<std::array<TH1D*, 3>>& pt_bin_hists,
                                     const char* canvas_name, const char* file_name) {
    const int n_pads = n_group_pt_bins + 1;
    const int n_columns = (n_pads + 1) / 2;
    TCanvas* canvas = new TCanvas(canvas_name, canvas_name, 600 * n_columns, 1000);
    canvas->Divide(n_columns, 2);
    canvas->cd(1);
    drawRegionHists(all_pt_hists, "all p_{T}");
    for (size_t i_pt_bin = 0; i_pt_bin < n_group_pt_bins; i_pt_bin++) {
      canvas->cd(i_pt_bin + 2);
      drawRegionHists(pt_bin_hists[i_pt_bin], ptBinLabel(i_pt_bin));
    }
    canvas->SaveAs(Form("%s/%s", output_dir.Data(), file_name));
  };
  drawPtBinnedRegionHists(hit_doublet_delta_phi_hists, pt_bin_doublet_delta_phi_hists,
                          "c_gen1_electron_hit_doublet_delta_phi", "gen1_electron_hit_doublet_delta_phi.pdf");
  drawPtBinnedRegionHists(hit_triplet_delta_cot_theta_hists, pt_bin_triplet_delta_cot_theta_hists,
                          "c_gen1_electron_hit_triplet_delta_cot_theta", "gen1_electron_hit_triplet_delta_cot_theta.pdf");

  // --- loopers: sum of the cot(theta) differences before and after the turn --
  TCanvas* looper_cot_theta_canvas = new TCanvas("c_gen1_electron_looper_delta_cot_theta",
                                                 "Conversion e+/e- looper delta cot theta", 800, 600);
  looper_sum_delta_cot_theta_hist->SetLineColor(kBlack);
  looper_delta_cot_theta_before_hist->SetLineColor(kBlue);
  looper_delta_cot_theta_before_hist->SetLineStyle(2);
  looper_sum_delta_cot_theta_hist->SetMinimum(0.5);
  looper_sum_delta_cot_theta_hist->SetMaximum(
      2 * std::max(looper_sum_delta_cot_theta_hist->GetMaximum(), looper_delta_cot_theta_before_hist->GetMaximum()));
  looper_sum_delta_cot_theta_hist->Draw("HIST");
  looper_delta_cot_theta_before_hist->Draw("HIST SAME");
  TLegend* looper_cot_theta_legend = new TLegend(0.58, 0.75, 0.88, 0.88);
  looper_cot_theta_legend->AddEntry(looper_sum_delta_cot_theta_hist, "#Deltacot#theta_{before} + #Deltacot#theta_{after}", "l");
  looper_cot_theta_legend->AddEntry(looper_delta_cot_theta_before_hist, "#Deltacot#theta_{before}", "l");
  looper_cot_theta_legend->Draw();
  gPad->SetLogy(1);  // narrow core, long tails
  looper_cot_theta_canvas->SaveAs(Form("%s/gen1_electron_looper_delta_cot_theta.pdf", output_dir.Data()));

  // --- loopers: pT estimate and circle centre r, after - before the turn ----
  TCanvas* looper_pt_canvas = new TCanvas("c_gen1_electron_looper_delta_pt", "Conversion e+/e- looper delta pT", 800, 600);
  looper_delta_pt_hist->SetLineColor(kBlue);
  looper_delta_pt_hist->Draw("HIST");
  gPad->SetLogy(1);
  looper_pt_canvas->SaveAs(Form("%s/gen1_electron_looper_delta_pt.pdf", output_dir.Data()));

  TCanvas* looper_center_canvas = new TCanvas("c_gen1_electron_looper_delta_center_r",
                                              "Conversion e+/e- looper delta circle centre r", 800, 600);
  looper_delta_center_r_hist->SetLineColor(kBlue);
  looper_delta_center_r_hist->Draw("HIST");
  gPad->SetLogy(1);
  looper_center_canvas->SaveAs(Form("%s/gen1_electron_looper_delta_center_r.pdf", output_dir.Data()));

  // --- production vertex r -------------------------------------------------
  TCanvas* r_canvas = new TCanvas("c_gen1_electron_vertex_r", "Conversion e+/e- vertex r", 800, 600);
  vertex_r_hist->SetLineColor(kBlue);
  vertex_r_hist->Draw("HIST");
  gPad->SetLogy(1);
  r_canvas->SaveAs(Form("%s/gen1_electron_vertex_r.pdf", output_dir.Data()));

  // --- production vertex r vs z ---------------------------------------------
  TCanvas* rz_canvas = new TCanvas("c_gen1_electron_vertex_rz", "Conversion e+/e- vertex r/z", 1000, 600);
  rz_canvas->SetRightMargin(0.14);
  vertex_rz_hist->GetZaxis()->SetTitleOffset(1.3);
  vertex_rz_hist->Draw("COLZ");
  gPad->SetLogz(1);
  rz_canvas->SaveAs(Form("%s/gen1_electron_vertex_rz.pdf", output_dir.Data()));
}
