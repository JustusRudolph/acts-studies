#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../../Utils/BasicFunctions.h"

constexpr double kSpeedOfLight = 299.792458;  // mm/ns

namespace {

// event, vertex_primary, vertex_secondary, particle, generation, sub_particle
using ParticleKey = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;

// a simulated hit; momentum in GeV, time in ns
struct SimHit {
  int index;  // position of the hit along the particle trajectory
  float x, y, z, t;
  float px, py, pz;
  double r() const { return std::hypot(x, y); }
  double phi() const { return std::atan2(y, x); }
  double pt() const { return std::hypot(px, py); }
  double phiDirection() const { return std::atan2(py, px); }
  bool isOutward() const { return x * px + y * py > 0; }
};

struct Electron {
  int pdg;
  std::vector<SimHit> layer_hits;  // hits in the studied layer, merged to one per crossing
};

// angle difference wrapped to [-pi, pi)
double deltaPhi(double phi_to, double phi_from) {
  return std::remainder(phi_to - phi_from, 2 * M_PI);
}

// hits of the same crossing (module overlaps) are less than merge_ds apart along the
// trajectory, a new crossing needs at least the short arc outside or inside the layer;
// keeps the first hit of each crossing
std::vector<SimHit> mergeCrossings(std::vector<SimHit> hits, double merge_ds) {
  std::sort(hits.begin(), hits.end(), [](const SimHit& a, const SimHit& b) { return a.index < b.index; });
  std::vector<SimHit> crossings;
  for (const SimHit& hit : hits) {
    if (!crossings.empty() && (hit.t - crossings.back().t) * kSpeedOfLight < merge_ds) continue;
    crossings.push_back(hit);
  }
  return crossings;
}

// k-th crossing: number of earlier crossings in the same direction (outward/inward),
// i.e. the turn relative to the first crossing in that direction
std::vector<int> turnIndices(const std::vector<SimHit>& crossings) {
  std::vector<int> turns;
  int n_out = 0, n_in = 0;
  for (const SimHit& hit : crossings) turns.push_back(hit.isOutward() ? n_out++ : n_in++);
  return turns;
}

// z spacing to the previous crossing in the same direction, two crossings back;
// NAN if that one is in the other direction (a crossing was missed)
double dzPerTurn(const std::vector<SimHit>& crossings, size_t k) {
  if (k < 2 || crossings[k].isOutward() != crossings[k - 2].isOutward()) return NAN;
  return crossings[k].z - crossings[k - 2].z;
}

// pitch of the helix [mm] from the momentum [GeV] in a field [T]: 2 pi R cot(theta)
double helixPitch(const SimHit& hit, double b_field) {
  return 2 * M_PI * hit.pz / (0.299792458e-3 * b_field);
}

double quantile(std::vector<double> values, double fraction) {
  if (values.empty()) return NAN;
  std::sort(values.begin(), values.end());
  const double position = fraction * (values.size() - 1);
  const size_t lower = size_t(position);
  if (lower + 1 >= values.size()) return values.back();
  return values[lower] + (position - lower) * (values[lower + 1] - values[lower]);
}

// symmetric plot range covering 99.5 % of |values|, with some margin
double plotRange(const std::vector<double>& values, double minimum) {
  std::vector<double> abs_values;
  for (double value : values) abs_values.push_back(std::abs(value));
  return std::max(minimum, 1.2 * quantile(abs_values, 0.995));
}

// 16 and 84 % quantile of values in each x bin, drawn as dashed lines
void drawQuantileLines(const std::vector<double>& x, const std::vector<double>& values, int x_min, int x_max) {
  TGraph* low = new TGraph();
  TGraph* high = new TGraph();
  for (int x_bin = x_min; x_bin <= x_max; x_bin++) {
    std::vector<double> in_bin;
    for (size_t i_value = 0; i_value < values.size(); i_value++) {
      if (int(x[i_value]) == x_bin) in_bin.push_back(values[i_value]);
    }
    if (in_bin.size() < 5) continue;
    low->AddPoint(x_bin, quantile(in_bin, 0.16));
    high->AddPoint(x_bin, quantile(in_bin, 0.84));
  }
  for (TGraph* graph : {low, high}) {
    graph->SetLineStyle(2);
    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);
    graph->Draw("L SAME");
  }
}

TH1D* fillHist(const char* name, const char* title, const std::vector<double>& values, double range, int n_bins) {
  TH1D* hist = new TH1D(name, title, n_bins, -range, range);
  for (double value : values) hist->Fill(value);
  return hist;
}

TH2D* fillHist2D(const char* name, const char* title, const std::vector<double>& x, const std::vector<double>& y,
                 int x_min, int x_max, double y_range, int n_bins_y) {
  TH2D* hist = new TH2D(name, title, x_max - x_min + 1, x_min - 0.5, x_max + 0.5, n_bins_y, -y_range, y_range);
  for (size_t i_value = 0; i_value < x.size(); i_value++) hist->Fill(x[i_value], y[i_value]);
  return hist;
}

void printQuantiles(const char* label, const std::vector<double>& values) {
  std::cout << Form("  %-48s 16/50/84 %%: %8.2f %8.2f %8.2f  (%zu entries)", label, quantile(values, 0.16),
                    quantile(values, 0.5), quantile(values, 0.84), values.size())
            << std::endl;
}

}  // namespace

// Real-particle version of the looper toy (scratch looper_toy.py / looper_turns.py):
// e+/e- that cross one barrel layer several times, by default the first outer
// tracker layer (volume 50, layer 3, r = 443 mm). Per crossing of that layer the
// true hit position, z, direction phi of the momentum, and whether the particle
// goes outward or inward; hits of the same crossing (module overlaps) are merged.
//
//   z spacing:   dz to the next crossing (alternates between the arc outside and the
//                arc inside the layer), and dz per turn, to the previous crossing in
//                the same direction; dz per turn is also compared to the pitch of the
//                helix 2 pi R cot(theta) from the true momentum at the earlier crossing
//   ddz:         consecutive dz(k+1) - dz(k), and same arc dz(k+2) - dz(k)
//   per turn:    r*delta phi of the hit and change of the direction phi, both relative
//                to the first crossing in the same direction
//
// Two figures in figures/looper_crossings/<input_dir>/, as for the toy; the top row of
// each shows one example particle, by default the one with the most crossings (event
// and barcode are printed for the event display).
//
// Usage (from pythia_studies):
//   root -l -b -q 'macros/looper_layer_crossings.C("gamma_conversion/10000ev_seed_1")'
void looper_layer_crossings(const TString input_dir="gamma_conversion/10000ev_seed_1",
                            const bool onSTBC=false,
                            const std::tuple<unsigned, unsigned> volume_layer={50, 3},  // OT layer, r = 443 mm
                            const unsigned min_crossings=4,
                            const unsigned max_turns=10,      // crossings beyond 2 * max_turns are ignored
                            const double merge_ds=10.,        // mm along the trajectory, same crossing
                            const double b_field=2.,          // T
                            const bool only_secondaries=false,  // generation > 0
                            const int example_event=-1,       // with example_particle; -1: most crossings
                            const unsigned example_particle=0) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetPaperSize(29, 17.7);  // cm, aspect ratio of the canvases; larger than A4 falls back to A4
  TH1::AddDirectory(false);
  const auto [layer_volume, layer_layer] = volume_layer;
  const unsigned max_crossings = 2 * max_turns;

  // --- e+/e- from the particles -----------------------------------------------
  const TString sim_filename = Utils::outputPathBase(onSTBC) + input_dir + "/particles_simulation.root";
  TFile* sim_file = TFile::Open(sim_filename);
  if (!sim_file || sim_file->IsZombie()) {
    std::cout << "Error: Could not open " << sim_filename << std::endl;
    return;
  }
  TTree* sim_tree = (TTree*)sim_file->Get("particles");
  if (!sim_tree) {
    std::cout << "Error: particles tree not found in " << sim_filename << std::endl;
    return;
  }
  unsigned event_id;
  std::vector<unsigned>* vertex_primary = nullptr;
  std::vector<unsigned>* vertex_secondary = nullptr;
  std::vector<unsigned>* particle = nullptr;
  std::vector<unsigned>* generation = nullptr;
  std::vector<unsigned>* sub_particle = nullptr;
  std::vector<int>* particle_type = nullptr;
  sim_tree->SetBranchStatus("*", false);
  for (const char* name : {"event_id", "vertex_primary", "vertex_secondary", "particle", "generation",
                           "sub_particle", "particle_type"}) {
    sim_tree->SetBranchStatus(name, true);
  }
  sim_tree->SetBranchAddress("event_id",         &event_id);
  sim_tree->SetBranchAddress("vertex_primary",   &vertex_primary);
  sim_tree->SetBranchAddress("vertex_secondary", &vertex_secondary);
  sim_tree->SetBranchAddress("particle",         &particle);
  sim_tree->SetBranchAddress("generation",       &generation);
  sim_tree->SetBranchAddress("sub_particle",     &sub_particle);
  sim_tree->SetBranchAddress("particle_type",    &particle_type);

  std::map<ParticleKey, Electron> electrons;
  const Long64_t n_events = sim_tree->GetEntries();
  for (Long64_t i_ev = 0; i_ev < n_events; i_ev++) {
    sim_tree->GetEntry(i_ev);
    for (size_t i_mcp = 0; i_mcp < particle_type->size(); i_mcp++) {
      if (std::abs(particle_type->at(i_mcp)) != 11) continue;
      if (only_secondaries && generation->at(i_mcp) == 0) continue;
      electrons[{event_id, vertex_primary->at(i_mcp), vertex_secondary->at(i_mcp), particle->at(i_mcp),
                 generation->at(i_mcp), sub_particle->at(i_mcp)}] = {particle_type->at(i_mcp), {}};
    }
  }
  sim_file->Close();
  delete sim_file;

  // --- hits: the studied layer for all e+/e-, all layers for the example -------
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
  float tx, ty, tz, tt, tpx, tpy, tpz;
  hits_tree->SetBranchStatus("*", false);
  for (const char* name : {"event_id", "barcode_vertex_primary", "barcode_vertex_secondary", "barcode_particle",
                           "barcode_generation", "barcode_sub_particle", "volume_id", "layer_id", "index",
                           "tx", "ty", "tz", "tt", "tpx", "tpy", "tpz"}) {
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
  hits_tree->SetBranchAddress("tt",                       &tt);
  hits_tree->SetBranchAddress("tpx",                      &tpx);
  hits_tree->SetBranchAddress("tpy",                      &tpy);
  hits_tree->SetBranchAddress("tpz",                      &tpz);
  const Long64_t n_hits = hits_tree->GetEntries();
  for (Long64_t i_hit = 0; i_hit < n_hits; i_hit++) {
    hits_tree->GetEntry(i_hit);
    if (volume_id != layer_volume || layer_id != layer_layer) continue;
    const auto electron = electrons.find(
        {hit_event_id, hit_vertex_primary, hit_vertex_secondary, hit_particle, hit_generation, hit_sub_particle});
    if (electron == electrons.end()) continue;
    electron->second.layer_hits.push_back({index, tx, ty, tz, tt, tpx, tpy, tpz});
  }

  // --- crossings per particle ---------------------------------------------------
  std::vector<std::pair<ParticleKey, std::vector<SimHit>>> loopers;  // crossings, at most max_crossings
  double layer_r_sum = 0;
  unsigned long n_layer_hits = 0;
  for (auto& [key, electron] : electrons) {
    if (electron.layer_hits.empty()) continue;
    for (const SimHit& hit : electron.layer_hits) layer_r_sum += hit.r();
    n_layer_hits += electron.layer_hits.size();
    std::vector<SimHit> crossings = mergeCrossings(electron.layer_hits, merge_ds);
    if (crossings.size() < min_crossings) continue;
    if (crossings.size() > max_crossings) crossings.resize(max_crossings);
    loopers.emplace_back(key, std::move(crossings));
  }
  if (loopers.empty()) {
    std::cout << "No e+/e- with >= " << min_crossings << " crossings of volume " << layer_volume << ", layer "
              << layer_layer << std::endl;
    return;
  }
  const double layer_r = layer_r_sum / n_layer_hits;

  // example particle: requested, else the one with the most crossings
  size_t i_example = 0;
  for (size_t i_looper = 0; i_looper < loopers.size(); i_looper++) {
    const ParticleKey& key = loopers[i_looper].first;
    if (example_event >= 0) {
      if (std::get<0>(key) == unsigned(example_event) && std::get<3>(key) == example_particle) i_example = i_looper;
    } else if (loopers[i_looper].second.size() > loopers[i_example].second.size()) {
      i_example = i_looper;
    }
  }
  const ParticleKey example_key = loopers[i_example].first;
  const std::vector<SimHit>& example = loopers[i_example].second;
  const int example_charge = electrons[example_key].pdg == 11 ? -1 : 1;

  // all hits of the example particle, for the transverse view
  std::vector<SimHit> example_all_hits;
  for (Long64_t i_hit = 0; i_hit < n_hits; i_hit++) {
    hits_tree->GetEntry(i_hit);
    if (ParticleKey{hit_event_id, hit_vertex_primary, hit_vertex_secondary, hit_particle, hit_generation,
                    hit_sub_particle} != example_key) continue;
    example_all_hits.push_back({index, tx, ty, tz, tt, tpx, tpy, tpz});
  }
  hits_file->Close();
  delete hits_file;

  // --- ensemble quantities ------------------------------------------------------
  std::vector<double> dz_turn, dz_turn_minus_pitch;
  std::vector<double> turn_of, rphi_change, direction_change;  // per crossing with turn > 0
  std::vector<double> ddz_consecutive, ddz_same_arc, k_consecutive, k_same_arc;
  std::map<unsigned, unsigned> n_crossings_count;
  for (const auto& [key, crossings] : loopers) {
    n_crossings_count[crossings.size()]++;
    const std::vector<int> turns = turnIndices(crossings);
    std::vector<double> dz;
    for (size_t k = 0; k < crossings.size(); k++) {
      if (k + 1 < crossings.size()) dz.push_back(crossings[k + 1].z - crossings[k].z);
      const double dz_k = dzPerTurn(crossings, k);
      if (std::isfinite(dz_k)) {
        dz_turn.push_back(dz_k);
        dz_turn_minus_pitch.push_back(dz_k - helixPitch(crossings[k - 2], b_field));
      }
      if (turns[k] == 0) continue;
      // first crossing in the same direction
      const SimHit& first = *std::find_if(crossings.begin(), crossings.end(), [&](const SimHit& hit) {
        return hit.isOutward() == crossings[k].isOutward();
      });
      turn_of.push_back(turns[k]);
      rphi_change.push_back(crossings[k].r() * deltaPhi(crossings[k].phi(), first.phi()));
      direction_change.push_back(deltaPhi(crossings[k].phiDirection(), first.phiDirection()) * 180 / M_PI);
    }
    for (size_t k = 0; k + 1 < dz.size(); k++) {
      ddz_consecutive.push_back(dz[k + 1] - dz[k]);
      k_consecutive.push_back(k);
      if (k + 2 < dz.size()) {
        ddz_same_arc.push_back(dz[k + 2] - dz[k]);
        k_same_arc.push_back(k);
      }
    }
  }

  std::cout << loopers.size() << " e+/e- with >= " << min_crossings << " crossings of volume " << layer_volume
            << ", layer " << layer_layer << " (mean r = " << Form("%.1f", layer_r) << " mm) in " << n_events
            << " events; crossings per particle (capped at " << max_crossings << "):";
  for (const auto& [n_crossings, count] : n_crossings_count) std::cout << " " << n_crossings << ": " << count << ",";
  std::cout << std::endl;
  std::cout << "Example: event " << std::get<0>(example_key) << ", barcode (" << std::get<1>(example_key) << ", "
            << std::get<2>(example_key) << ", " << std::get<3>(example_key) << ", " << std::get<4>(example_key) << ", "
            << std::get<5>(example_key) << "), pdg " << electrons[example_key].pdg << ", " << example.size()
            << " crossings, pT at the first " << Form("%.1f", 1e3 * example.front().pt()) << " MeV, cot theta "
            << Form("%.3f", example.front().pz / example.front().pt()) << std::endl;
  printQuantiles("dz per turn [mm]", dz_turn);
  printQuantiles("dz per turn - pitch from momentum [mm]", dz_turn_minus_pitch);
  printQuantiles("ddz consecutive [mm]", ddz_consecutive);
  printQuantiles("ddz same arc [mm]", ddz_same_arc);
  for (int turn = 1; turn < int(max_turns); turn++) {
    std::vector<double> rphi_turn, direction_turn;
    for (size_t i_value = 0; i_value < turn_of.size(); i_value++) {
      if (int(turn_of[i_value]) != turn) continue;
      rphi_turn.push_back(rphi_change[i_value]);
      direction_turn.push_back(direction_change[i_value]);
    }
    if (rphi_turn.empty()) continue;
    printQuantiles(Form("turn %d: r*delta phi of the hit [mm]", turn), rphi_turn);
    printQuantiles(Form("turn %d: direction phi change [deg]", turn), direction_turn);
  }

  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/pythia_studies/figures/looper_crossings/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);
  const TString layer_tag = Form("vol%u_layer%u", layer_volume, layer_layer);
  const TString layer_text = Form("volume %u, layer %u (r = %.0f mm)", layer_volume, layer_layer, layer_r);
  const std::vector<int> example_turns = turnIndices(example);
  auto turnColor = [&](int turn) { return TColor::GetColorPalette(turn * (gStyle->GetNumberOfColors() - 1) / std::max(1u, max_turns - 1)); };

  // --- figure 1: example particle and per-turn spread ---------------------------
  TCanvas* main_canvas = new TCanvas("c_looper_crossings", "Looper layer crossings", 1800, 1100);
  main_canvas->Divide(3, 2);

  // transverse view: all hits, the layer, the circle from the momentum at the first crossing
  main_canvas->cd(1);
  gPad->SetLeftMargin(0.14);
  double x_min = 1e9, x_max = -1e9, y_min = 1e9, y_max = -1e9;
  for (const SimHit& hit : example_all_hits) {
    x_min = std::min<double>(x_min, hit.x); x_max = std::max<double>(x_max, hit.x);
    y_min = std::min<double>(y_min, hit.y); y_max = std::max<double>(y_max, hit.y);
  }
  const double half_size = 0.55 * std::max(x_max - x_min, y_max - y_min) + 20;
  const double x_mid = 0.5 * (x_min + x_max), y_mid = 0.5 * (y_min + y_max);
  TH2D* xy_frame = new TH2D("xy_frame", Form("example %s, transverse view;x (mm);y (mm)", example_charge < 0 ? "e^{-}" : "e^{+}"),
                            10, x_mid - half_size, x_mid + half_size, 10, y_mid - half_size, y_mid + half_size);
  xy_frame->Draw();
  TEllipse* layer_circle = new TEllipse(0, 0, layer_r, layer_r);
  layer_circle->SetFillStyle(0);
  layer_circle->Draw();
  const SimHit& first = example.front();
  const double radius = first.pt() / (0.299792458e-3 * b_field);
  // centre: v x B points along (py, -px) for a positive charge
  TEllipse* helix_circle = new TEllipse(first.x + example_charge * radius * first.py / first.pt(),
                                        first.y - example_charge * radius * first.px / first.pt(), radius, radius);
  helix_circle->SetFillStyle(0);
  helix_circle->SetLineStyle(2);
  helix_circle->SetLineColor(kRed);
  helix_circle->Draw();
  TGraph* all_hits_graph = new TGraph();
  for (const SimHit& hit : example_all_hits) all_hits_graph->AddPoint(hit.x, hit.y);
  all_hits_graph->SetMarkerStyle(20);
  all_hits_graph->SetMarkerSize(0.4);
  all_hits_graph->SetMarkerColor(kGray + 1);
  all_hits_graph->Draw("P SAME");
  for (size_t k = 0; k < example.size(); k++) {
    TGraph* point = new TGraph();
    point->AddPoint(example[k].x, example[k].y);
    point->SetMarkerStyle(example[k].isOutward() ? 20 : 21);
    point->SetMarkerColor(turnColor(example_turns[k]));
    point->Draw("P SAME");
  }
  TLegend* xy_legend = new TLegend(0.15, 0.80, 0.60, 0.89);
  xy_legend->SetTextSize(0.03);
  xy_legend->AddEntry(helix_circle, "circle from p at the first crossing", "l");
  xy_legend->AddEntry(all_hits_graph, "all hits", "p");
  xy_legend->Draw();

  // hits on the layer: r*phi vs z, colour = turn, circle = outward, square = inward
  main_canvas->cd(2);
  gPad->SetLeftMargin(0.14);
  TGraph* rphi_z_frame = new TGraph();
  // r*phi relative to the first crossing, so that hits around phi = +-pi stay together
  auto exampleRPhi = [&](const SimHit& hit) { return hit.r() * deltaPhi(hit.phi(), example.front().phi()); };
  for (const SimHit& hit : example) rphi_z_frame->AddPoint(hit.z, exampleRPhi(hit));
  rphi_z_frame->SetTitle("example: hits on the layer (#bullet out, #blacksquare in, colour = turn);z (mm);"
                         "r#Delta#phi to the first crossing (mm)");
  rphi_z_frame->SetMarkerStyle(1);
  rphi_z_frame->Draw("AP");
  for (size_t k = 0; k < example.size(); k++) {
    TGraph* point = new TGraph();
    point->AddPoint(example[k].z, exampleRPhi(example[k]));
    point->SetMarkerStyle(example[k].isOutward() ? 20 : 21);
    point->SetMarkerSize(1.2);
    point->SetMarkerColor(turnColor(example_turns[k]));
    point->Draw("P SAME");
  }

  // z spacing of the example
  main_canvas->cd(3);
  gPad->SetLeftMargin(0.14);
  TGraph* dz_next_graph = new TGraph();
  TGraph* dz_turn_graph = new TGraph();
  for (size_t k = 1; k < example.size(); k++) {
    dz_next_graph->AddPoint(k, example[k].z - example[k - 1].z);
    const double dz_k = dzPerTurn(example, k);
    if (std::isfinite(dz_k)) dz_turn_graph->AddPoint(k, dz_k);
  }
  dz_next_graph->SetTitle("example: z spacing;crossing;dz (mm)");
  dz_next_graph->SetMarkerStyle(20);
  dz_next_graph->SetLineColor(kBlue);
  dz_next_graph->SetMarkerColor(kBlue);
  dz_turn_graph->SetMarkerStyle(21);
  dz_turn_graph->SetMarkerColor(kOrange + 7);
  dz_next_graph->Draw("APL");
  dz_turn_graph->Draw("P SAME");
  TLegend* dz_legend = new TLegend(0.50, 0.80, 0.89, 0.89);
  dz_legend->SetTextSize(0.03);
  dz_legend->AddEntry(dz_next_graph, "to the next crossing", "lp");
  dz_legend->AddEntry(dz_turn_graph, "per turn (same direction)", "p");
  dz_legend->Draw();

  // ensemble: dz per turn, and its difference to the pitch from the momentum
  main_canvas->cd(4);
  const double dz_range = plotRange(dz_turn, 10);
  TH1D* dz_turn_hist = fillHist("dz_turn_hist", Form("%zu e^{#pm}: dz per turn;dz (mm);Entries", loopers.size()),
                                dz_turn, dz_range, 150);
  TH1D* dz_pitch_hist = fillHist("dz_pitch_hist", "", dz_turn_minus_pitch, dz_range, 150);
  dz_turn_hist->SetLineColor(kBlue);
  dz_pitch_hist->SetLineColor(kRed);
  dz_turn_hist->SetMinimum(0.5);
  dz_turn_hist->SetMaximum(2 * std::max(dz_turn_hist->GetMaximum(), dz_pitch_hist->GetMaximum()));
  dz_turn_hist->Draw("HIST");
  dz_pitch_hist->Draw("HIST SAME");
  gPad->SetLogy(1);
  TLegend* pitch_legend = new TLegend(0.45, 0.78, 0.89, 0.89);
  pitch_legend->SetTextSize(0.03);
  pitch_legend->AddEntry(dz_turn_hist, "dz per turn", "l");
  pitch_legend->AddEntry(dz_pitch_hist, "dz per turn - 2#piRcot#theta (p at the earlier crossing)", "l");
  pitch_legend->Draw();

  // ensemble: hit position and direction relative to the first crossing in the same direction
  for (int i_pad : {5, 6}) {
    main_canvas->cd(i_pad);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.14);
    const bool isPosition = i_pad == 5;
    const std::vector<double>& values = isPosition ? rphi_change : direction_change;
    TH2D* turn_hist = fillHist2D(isPosition ? "rphi_turn_hist" : "direction_turn_hist",
                                 isPosition ? "hit position along the layer (same-direction crossings);turn;"
                                              "r#Delta#phi of the hit vs first turn (mm)"
                                            : "direction #phi at the crossing (same-direction crossings);turn;"
                                              "direction #phi change vs first turn (deg)",
                                 turn_of, values, 1, max_turns - 1, plotRange(values, 1), 100);
    turn_hist->SetMinimum(1);
    turn_hist->Draw("COLZ");
    gPad->SetLogz(1);
    drawQuantileLines(turn_of, values, 1, max_turns - 1);
  }
  main_canvas->cd(0);
  main_canvas->SaveAs(Form("%s/looper_crossings_%s.pdf", output_dir.Data(), layer_tag.Data()));

  // --- figure 2: ddz -------------------------------------------------------------
  TCanvas* ddz_canvas = new TCanvas("c_looper_crossings_ddz", "Looper layer crossings ddz", 1800, 1100);
  ddz_canvas->Divide(3, 2);

  ddz_canvas->cd(1);
  gPad->SetLeftMargin(0.14);
  TGraph* example_consecutive = new TGraph();
  TGraph* example_same_arc = new TGraph();
  for (size_t k = 0; k + 2 < example.size(); k++) {
    const double dz_k = example[k + 1].z - example[k].z;
    example_consecutive->AddPoint(k + 1, (example[k + 2].z - example[k + 1].z) - dz_k);
    if (k + 3 < example.size()) example_same_arc->AddPoint(k + 2, (example[k + 3].z - example[k + 2].z) - dz_k);
  }
  example_consecutive->SetTitle(Form("example (%zu crossings);k;ddz (mm)", example.size()));
  example_consecutive->SetMarkerStyle(20);
  example_consecutive->SetMarkerColor(kBlue);
  example_consecutive->SetLineColor(kBlue);
  example_same_arc->SetMarkerStyle(21);
  example_same_arc->SetMarkerColor(kOrange + 7);
  example_consecutive->Draw("APL");
  example_same_arc->Draw("P SAME");
  TLine* zero_line = new TLine(example_consecutive->GetXaxis()->GetXmin(), 0,
                               example_consecutive->GetXaxis()->GetXmax(), 0);
  zero_line->Draw();
  TLegend* ddz_legend = new TLegend(0.50, 0.80, 0.89, 0.89);
  ddz_legend->SetTextSize(0.03);
  ddz_legend->AddEntry(example_consecutive, "consecutive: dz(k+1) - dz(k)", "lp");
  ddz_legend->AddEntry(example_same_arc, "same arc: dz(k+2) - dz(k)", "p");
  ddz_legend->Draw();

  const double ddz_range = plotRange(ddz_consecutive, 5);
  for (int i_pad : {2, 3}) {
    ddz_canvas->cd(i_pad);
    const bool isConsecutive = i_pad == 2;
    const std::vector<double>& values = isConsecutive ? ddz_consecutive : ddz_same_arc;
    TH1D* hist = fillHist(isConsecutive ? "ddz_consecutive_hist" : "ddz_same_arc_hist",
                          Form("%s: median %.2f, 16-84 %%: [%.2f, %.2f] mm;ddz (mm);Entries",
                               isConsecutive ? "consecutive dz(k+1) - dz(k)" : "same arc dz(k+2) - dz(k)",
                               quantile(values, 0.5), quantile(values, 0.16), quantile(values, 0.84)),
                          values, ddz_range, 150);
    hist->SetLineColor(kBlue);
    hist->SetMinimum(0.5);
    hist->Draw("HIST");
    gPad->SetLogy(1);
  }

  // half 16-84 % width of ddz per k
  ddz_canvas->cd(4);
  gPad->SetLeftMargin(0.14);
  TGraph* width_graphs[2] = {new TGraph(), new TGraph()};
  for (int i_kind = 0; i_kind < 2; i_kind++) {
    const std::vector<double>& ks = i_kind == 0 ? k_consecutive : k_same_arc;
    const std::vector<double>& values = i_kind == 0 ? ddz_consecutive : ddz_same_arc;
    for (int k = 0; k + 2 < int(max_crossings); k++) {
      std::vector<double> in_bin;
      for (size_t i_value = 0; i_value < values.size(); i_value++) {
        if (int(ks[i_value]) == k) in_bin.push_back(values[i_value]);
      }
      if (in_bin.size() < 5) continue;
      width_graphs[i_kind]->AddPoint(k, 0.5 * (quantile(in_bin, 0.84) - quantile(in_bin, 0.16)));
    }
    width_graphs[i_kind]->SetMarkerStyle(20 + i_kind);
    width_graphs[i_kind]->SetMarkerColor(i_kind == 0 ? kBlue : kOrange + 7);
    width_graphs[i_kind]->SetLineColor(i_kind == 0 ? kBlue : kOrange + 7);
  }
  width_graphs[0]->SetTitle("spread of ddz per k;k;half 16-84 % width of ddz (mm)");
  width_graphs[0]->SetMinimum(0);
  width_graphs[0]->SetMaximum(1.2 * std::max(TMath::MaxElement(width_graphs[0]->GetN(), width_graphs[0]->GetY()),
                                             width_graphs[1]->GetN() ? TMath::MaxElement(width_graphs[1]->GetN(), width_graphs[1]->GetY()) : 0.));
  width_graphs[0]->Draw("APL");
  width_graphs[1]->Draw("PL SAME");
  TLegend* width_legend = new TLegend(0.60, 0.15, 0.89, 0.27);
  width_legend->AddEntry(width_graphs[0], "consecutive", "lp");
  width_legend->AddEntry(width_graphs[1], "same arc", "lp");
  width_legend->Draw();

  for (int i_pad : {5, 6}) {
    ddz_canvas->cd(i_pad);
    gPad->SetRightMargin(0.14);
    const bool isConsecutive = i_pad == 5;
    TH2D* hist = fillHist2D(isConsecutive ? "ddz_consecutive_k_hist" : "ddz_same_arc_k_hist",
                            isConsecutive ? "ddz vs k, consecutive;k;ddz (mm)" : "ddz vs k, same arc;k;ddz (mm)",
                            isConsecutive ? k_consecutive : k_same_arc, isConsecutive ? ddz_consecutive : ddz_same_arc,
                            0, max_crossings - 3, ddz_range, 120);
    hist->SetMinimum(1);
    hist->Draw("COLZ");
    gPad->SetLogz(1);
  }
  ddz_canvas->cd(0);
  ddz_canvas->SaveAs(Form("%s/looper_crossings_%s_ddz.pdf", output_dir.Data(), layer_tag.Data()));
  std::cout << "Layer: " << layer_text << std::endl;
}
