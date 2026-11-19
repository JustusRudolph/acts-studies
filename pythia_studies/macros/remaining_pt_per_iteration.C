// pT distribution of the particles that are still not reconstructed after each
// tracking iteration, i.e. not found by any track from that iteration or an
// earlier one.
//
// Left: all simulated particles with at least one hit.
// Right: the selected particles (particles_digitized_selected.root).
//
// A particle counts as found when it is the majority particle of a track and
// owns at least `matching_ratio` of that track's measurements, the ACTS truth
// matching criterion. The per-iteration track summaries only name the majority
// particle when it is a selected one, so the hit-level truth is taken from the
// track states of the merged, ambiguity-resolved tracks instead. Each merged
// track is attributed to its iteration by matching (event, nMeasurements,
// chi2Sum, fitted q/p) against the per-iteration ambiguity-resolved track
// summaries. Tracks that the final ambiguity resolution removed are not
// counted.
//
// Usage (from pythia_studies):
//   root -l -b -q 'macros/remaining_pt_per_iteration.C+("<run dir>", 5)'
// <run dir> holds particles_digitized_selected.root, the merged
// trackstates/tracksummary files and iter_<i>/. The simulated particles are
// read from <sim dir>/particles_simulation.root, where <sim dir> is the run dir
// itself for a full-chain run, or, for a reco_only run, the sibling directory
// named by --simdir in configs/script_command_line_args_used.log. Pass sim_dir
// to override.

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// dictionary for the per-state particle ids of the track states (needs ACLiC, i.e. the + when running)
#ifdef __ROOTCLING__
#pragma link C++ class std::vector<std::vector<unsigned int>>+;
#endif

namespace {

// event, vertex_primary, vertex_secondary, particle, generation, sub_particle
using ParticleKey = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;
// event, nMeasurements, chi2Sum bits, fitted q/p bits
using TrackKey = std::tuple<unsigned, unsigned, unsigned, unsigned>;

constexpr int kNotFound = 1 << 30;

struct ParticleInfo {
  float pt;
  bool isSelected = false;
  int foundInIteration = kNotFound;
};

unsigned floatBits(float value) {
  unsigned bits;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

// iteration 0 writes its files without the iter_0_ prefix
TString ambiSummaryFile(const TString& data_dir, int i_iter) {
  TString prefix = i_iter == 0 ? "" : Form("iter_%d_", i_iter);
  return Form("%s/iter_%d/tracksummary_%sambi_tracks_writer.root", data_dir.Data(), i_iter, prefix.Data());
}

// directory holding particles_simulation.root, see the header
TString findSimDir(const TString& data_dir) {
  if (!gSystem->AccessPathName(data_dir + "/particles_simulation.root")) return data_dir;
  std::ifstream log((data_dir + "/configs/script_command_line_args_used.log").Data());
  std::string token;
  while (log >> token) {
    if (token != "--simdir" || !(log >> token)) continue;
    TString sim_dir = TString(token).Strip(TString::kTrailing, '/');
    return TString(gSystem->GetDirName(data_dir)) + "/" + gSystem->BaseName(sim_dir.Data());
  }
  return "";
}

// (event, nMeasurements, chi2Sum, q/p) of every track in a track summary file,
// in the order of the file; returns false if the file cannot be read
bool readTrackKeys(const TString& path, std::vector<TrackKey>& keys, std::vector<unsigned>* track_nrs = nullptr) {
  TFile* file = TFile::Open(path);
  if (!file || file->IsZombie()) return false;
  TTree* tree = (TTree*) file->Get("tracksummary");
  tree->SetBranchStatus("*", false);
  unsigned event_nr;
  std::vector<unsigned>* track_nr = nullptr;
  std::vector<unsigned>* nMeasurements = nullptr;
  std::vector<float>* chi2Sum = nullptr;
  std::vector<float>* qop = nullptr;
  for (const char* name : {"event_nr", "track_nr", "nMeasurements", "chi2Sum", "eQOP_fit"}) {
    tree->SetBranchStatus(name, true);
  }
  tree->SetBranchAddress("event_nr", &event_nr);
  tree->SetBranchAddress("track_nr", &track_nr);
  tree->SetBranchAddress("nMeasurements", &nMeasurements);
  tree->SetBranchAddress("chi2Sum", &chi2Sum);
  tree->SetBranchAddress("eQOP_fit", &qop);
  for (Long64_t i_event = 0; i_event < tree->GetEntries(); i_event++) {
    tree->GetEntry(i_event);
    for (size_t i_track = 0; i_track < nMeasurements->size(); i_track++) {
      keys.emplace_back(event_nr, nMeasurements->at(i_track), floatBits(chi2Sum->at(i_track)),
                        floatBits(qop->at(i_track)));
      if (track_nrs) track_nrs->push_back(track_nr->at(i_track));
    }
  }
  file->Close();
  delete file;
  return true;
}

// logarithmic binning between pt_min and pt_max [GeV]
std::vector<double> logBins(int n_bins, double pt_min, double pt_max) {
  std::vector<double> edges(n_bins + 1);
  for (int i_bin = 0; i_bin <= n_bins; i_bin++) {
    edges[i_bin] = pt_min * std::pow(pt_max / pt_min, double(i_bin) / n_bins);
  }
  return edges;
}

}  // namespace

void remaining_pt_per_iteration(const std::string data_dir_arg,
                                const int n_iterations = 4,
                                const std::string sim_dir_arg = "",
                                const double matching_ratio = 0.5) {
  gStyle->SetOptStat(0);
  const TString data_dir = TString(data_dir_arg).Strip(TString::kTrailing, '/');
  const TString sim_dir = sim_dir_arg.empty() ? findSimDir(data_dir) : TString(sim_dir_arg);
  if (sim_dir.IsNull()) {
    std::cout << "Error: could not determine the simulation directory, pass sim_dir" << std::endl;
    return;
  }

  // --- simulated particles with hits, and which of them are selected -------
  std::map<ParticleKey, ParticleInfo> particles;
  Long64_t n_events = 0;
  auto readParticles = [&](const TString& path, bool isSelection) {
    TFile* file = TFile::Open(path);
    if (!file || file->IsZombie()) {
      std::cout << "Error: could not open " << path << std::endl;
      return false;
    }
    TTree* tree = (TTree*) file->Get("particles");
    tree->SetBranchStatus("*", false);
    unsigned event_id;
    std::vector<unsigned>* vertex_primary = nullptr;
    std::vector<unsigned>* vertex_secondary = nullptr;
    std::vector<unsigned>* particle = nullptr;
    std::vector<unsigned>* generation = nullptr;
    std::vector<unsigned>* sub_particle = nullptr;
    std::vector<float>* pt = nullptr;
    std::vector<int>* number_of_hits = nullptr;
    for (const char* name : {"event_id", "vertex_primary", "vertex_secondary", "particle", "generation",
                             "sub_particle", "pt", "number_of_hits"}) {
      tree->SetBranchStatus(name, true);
    }
    tree->SetBranchAddress("event_id", &event_id);
    tree->SetBranchAddress("vertex_primary", &vertex_primary);
    tree->SetBranchAddress("vertex_secondary", &vertex_secondary);
    tree->SetBranchAddress("particle", &particle);
    tree->SetBranchAddress("generation", &generation);
    tree->SetBranchAddress("sub_particle", &sub_particle);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("number_of_hits", &number_of_hits);

    long n_not_simulated = 0;
    for (Long64_t i_event = 0; i_event < tree->GetEntries(); i_event++) {
      tree->GetEntry(i_event);
      for (size_t i_part = 0; i_part < pt->size(); i_part++) {
        if (!isSelection && number_of_hits->at(i_part) < 1) continue;
        ParticleKey key{event_id, vertex_primary->at(i_part), vertex_secondary->at(i_part),
                        particle->at(i_part), generation->at(i_part), sub_particle->at(i_part)};
        if (!isSelection) {
          particles[key].pt = pt->at(i_part);
          continue;
        }
        auto it = particles.find(key);
        if (it == particles.end()) {
          n_not_simulated++;
          continue;
        }
        it->second.isSelected = true;
      }
    }
    if (!isSelection) n_events = tree->GetEntries();
    if (n_not_simulated) printf("Warning: %ld selected particles have no simulated particle with hits\n", n_not_simulated);
    file->Close();
    delete file;
    return true;
  };
  if (!readParticles(sim_dir + "/particles_simulation.root", false)) return;
  if (!readParticles(data_dir + "/particles_digitized_selected.root", true)) return;
  long n_selected = 0;
  for (const auto& [key, info] : particles) n_selected += info.isSelected;
  printf("Loaded %zu simulated particles with hits (%ld selected) from %lld events in %s\n",
         particles.size(), n_selected, n_events, sim_dir.Data());

  // --- iteration of each merged track ----------------------------------------
  std::map<TrackKey, int> iterationOfTrack;
  int n_read = 0;
  for (int i_iter = 0; i_iter < n_iterations; i_iter++) {
    std::vector<TrackKey> keys;
    if (!readTrackKeys(ambiSummaryFile(data_dir, i_iter), keys)) {
      std::cout << "Warning: could not open " << ambiSummaryFile(data_dir, i_iter) << ", stopping at "
                << i_iter << " iterations" << std::endl;
      break;
    }
    for (const TrackKey& key : keys) iterationOfTrack[key] = i_iter;
    n_read++;
  }
  std::vector<TrackKey> merged_keys;
  std::vector<unsigned> merged_track_nrs;
  if (!readTrackKeys(data_dir + "/tracksummary_ambi-tracks-merged.root", merged_keys, &merged_track_nrs)) {
    std::cout << "Error: could not open the merged track summary in " << data_dir << std::endl;
    return;
  }
  // (event, track_nr) -> iteration
  std::map<std::pair<unsigned, unsigned>, int> iterationOfMerged;
  long n_unattributed = 0;
  for (size_t i_track = 0; i_track < merged_keys.size(); i_track++) {
    auto it = iterationOfTrack.find(merged_keys[i_track]);
    if (it == iterationOfTrack.end()) {
      n_unattributed++;
      continue;
    }
    iterationOfMerged[{std::get<0>(merged_keys[i_track]), merged_track_nrs[i_track]}] = it->second;
  }
  printf("Attributed %zu of %zu merged tracks to an iteration\n", iterationOfMerged.size(), merged_keys.size());

  // --- majority particle of each merged track from its track states ---------
  TFile* states_file = TFile::Open(data_dir + "/trackstates_ambi-tracks-merged.root");
  if (!states_file || states_file->IsZombie()) {
    std::cout << "Error: could not open the merged track states in " << data_dir << std::endl;
    return;
  }
  TTree* states = (TTree*) states_file->Get("trackstates");
  states->SetBranchStatus("*", false);
  unsigned event_nr, track_nr;
  std::vector<int>* stateType = nullptr;
  std::vector<std::vector<unsigned>>* ids[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  const char* id_names[5] = {"particle_ids_vertex_primary", "particle_ids_vertex_secondary",
                             "particle_ids_particle", "particle_ids_generation", "particle_ids_sub_particle"};
  for (const char* name : {"event_nr", "track_nr", "stateType"}) states->SetBranchStatus(name, true);
  states->SetBranchAddress("event_nr", &event_nr);
  states->SetBranchAddress("track_nr", &track_nr);
  states->SetBranchAddress("stateType", &stateType);
  for (int i_id = 0; i_id < 5; i_id++) {
    states->SetBranchStatus(id_names[i_id], true);
    states->SetBranchAddress(id_names[i_id], &ids[i_id]);
  }

  std::vector<long> n_found(n_read, 0), n_found_selected(n_read, 0);
  long n_found_unknown = 0;
  std::map<ParticleKey, unsigned> hitCounts;
  for (Long64_t i_track = 0; i_track < states->GetEntries(); i_track++) {
    states->GetEntry(i_track);
    auto it_iter = iterationOfMerged.find({event_nr, track_nr});
    if (it_iter == iterationOfMerged.end()) continue;

    // every particle of every measurement counts, as in the ACTS truth matching
    hitCounts.clear();
    unsigned nMeasurements = 0;
    for (size_t i_state = 0; i_state < stateType->size(); i_state++) {
      if (stateType->at(i_state) != 0) continue;  // measurements only
      nMeasurements++;
      for (size_t i_part = 0; i_part < ids[2]->at(i_state).size(); i_part++) {
        hitCounts[{event_nr, ids[0]->at(i_state)[i_part], ids[1]->at(i_state)[i_part], ids[2]->at(i_state)[i_part],
                   ids[3]->at(i_state)[i_part], ids[4]->at(i_state)[i_part]}]++;
      }
    }
    const ParticleKey* majority = nullptr;
    unsigned nMajority = 0;
    for (const auto& [key, count] : hitCounts) {
      if (count > nMajority) {
        majority = &key;
        nMajority = count;
      }
    }
    if (!majority || nMajority < matching_ratio * nMeasurements) continue;  // fake

    auto it = particles.find(*majority);
    if (it == particles.end()) {
      n_found_unknown++;
      continue;
    }
    if (it->second.foundInIteration <= it_iter->second) continue;  // found before
    if (it->second.foundInIteration != kNotFound) {  // moved to an earlier iteration
      n_found[it->second.foundInIteration]--;
      n_found_selected[it->second.foundInIteration] -= it->second.isSelected;
    }
    it->second.foundInIteration = it_iter->second;
    n_found[it_iter->second]++;
    n_found_selected[it_iter->second] += it->second.isSelected;
  }
  states_file->Close();
  delete states_file;
  if (n_found_unknown) printf("Warning: %ld tracks have a majority particle without simulated hits\n", n_found_unknown);

  // --- histograms: before any iteration, then after each one ---------------
  std::vector<double> edges = logBins(60, 0.01, 20.);
  std::vector<TH1D*> h_all, h_selected;
  for (int i_hist = 0; i_hist <= n_read; i_hist++) {
    TString tag = i_hist == 0 ? TString("initial") : Form("iter_%d", i_hist - 1);
    h_all.push_back(new TH1D("h_all_" + tag, "", edges.size() - 1, edges.data()));
    h_selected.push_back(new TH1D("h_selected_" + tag, "", edges.size() - 1, edges.data()));
    for (const auto& [key, info] : particles) {
      if (info.foundInIteration < i_hist) continue;  // found in iteration i_hist - 1 or earlier
      h_all.back()->Fill(info.pt);
      if (info.isSelected) h_selected.back()->Fill(info.pt);
    }
    if (i_hist > 0) {
      printf("Iteration %d: %ld newly found (%ld selected), %.0f remaining (%.0f selected)\n", i_hist - 1,
             n_found[i_hist - 1], n_found_selected[i_hist - 1], h_all.back()->GetEntries(),
             h_selected.back()->GetEntries());
    }
  }

  // --- draw -----------------------------------------------------------------
  const std::vector<int> colors = {kBlack, kAzure + 1, kGreen + 2, kOrange + 7,
                                   kRed + 1, kMagenta + 1, kViolet + 2, kCyan + 2, kGray + 2};
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1600, 700);
  canvas->Divide(2, 1);
  TLegend* legend = new TLegend(0.55, 0.66, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  const TString titles[2] = {"All particles with hits", "Selected particles"};
  for (int i_pad = 0; i_pad < 2; i_pad++) {
    std::vector<TH1D*>& hists = i_pad == 0 ? h_all : h_selected;
    canvas->cd(i_pad + 1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.13);

    double y_max = 0.;
    for (TH1D* h : hists) {
      h->Scale(1. / n_events);
      y_max = std::max(y_max, h->GetMaximum());
    }
    for (size_t i_hist = 0; i_hist < hists.size(); i_hist++) {
      TH1D* h = hists[i_hist];
      h->SetLineColor(colors[i_hist % colors.size()]);
      h->SetLineWidth(2);
      h->SetTitle(titles[i_pad] + ";#it{p}_{T} [GeV/#it{c}];Not yet reconstructed particles / event");
      h->GetYaxis()->SetRangeUser(5e-3, 50 * y_max);  // headroom for the legend
      h->Draw(i_hist == 0 ? "HIST" : "HIST SAME");
      if (i_pad == 0) {
        legend->AddEntry(h, i_hist == 0 ? TString("Before iteration 0")
                                        : Form("After iteration %zu", i_hist - 1), "l");
      }
    }
    legend->Draw();
  }

  // output name from the run directory's name
  TString run_name = gSystem->BaseName(data_dir.Data());
  gSystem->mkdir("figures/remaining_pt_per_iteration", true);
  gSystem->mkdir("hists/remaining_pt_per_iteration", true);
  TString output_name = Form("%s_%d_iter", run_name.Data(), n_read);
  canvas->SaveAs(Form("figures/remaining_pt_per_iteration/%s.pdf", output_name.Data()));

  TFile* out_file = TFile::Open(Form("hists/remaining_pt_per_iteration/%s.root", output_name.Data()), "RECREATE");
  for (TH1D* h : h_all) h->Write();
  for (TH1D* h : h_selected) h->Write();
  out_file->Close();
  delete out_file;
  delete canvas;
}
