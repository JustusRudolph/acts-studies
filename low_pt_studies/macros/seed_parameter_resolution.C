#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

#include "../../Utils/BasicFunctions.h"

// Resolution of the seed parameters against the truth, per true pT, compared with the
// sigmas their covariance assigns. Which parameters is set by `pass`:
//  0: the estimated seed parameters, i.e. the CKF start
//     (<iter>_estimatedparameters.root, ACTS RootTrackParameterWriter, written after
//     the TOF extension / seed flip). The truth there is the particle at the first
//     measurement of the seed.
//  1: the CKF state after the seed: filtered on the last seed hit (in the CKF
//     direction; the TOF hit of an unflipped TOF extended seed), i.e. the seed fitted
//     through its own hits (<iter>_seedexitparameters.root, RootSeedExitWriter).
//  2: the start of the CKF second pass: the smoothed first pass state on the first
//     measurement of the first pass track (same file).
// For 1 and 2 the truth is the seed's majority particle at the hit the state sits on;
// seeds without that state (status != 0, e.g. deduplicated, or the last seed hit not
// taken by the CKF) or whose state sits on a hit of another particle are excluded.
//
// Seeds used:
//  - truthMatched (all seed hits with a truth particle belong to one particle; hits
//    without truth, i.e. noise, are ignored by that flag, and only those seeds have
//    truth parameters), and
//  - exactly own_hits of the three space points of the seed (seeds.root, same order as
//    the parameters) from that particle, from the truth of the measurements
//    (measurements_with_noise.root): 3 = a correct seed, 2 = two hits of the particle and
//    one noise hit (a hit of another particle leaves the seed without truth). A TOF hit
//    added by the TOF extension is not part of the triplet and not checked, and
//  - pass 0: parameters and truth on the same surface (|res_loc0|, |res_loc1| < 0.1 mm;
//    the truth is taken at the first measurement of the seed, which for a flipped seed
//    is not the one its parameters are bound to); pass 1, 2: a state on a hit of the
//    particle (status 0, t_on_particle); pass 1 also at the end of the seed (on its
//    third space point or TOF hit, with all three tracker hits taken), and
//  - the same charge as the truth (seeds running backwards in time are excluded).
// The numbers of excluded seeds per reason are in the last panel.
//
// Residual = estimate - truth, signed (loc0, loc1 in um; phi, theta in degrees; q/p
// relative: (|q/p| - |q/p_true|) / |q/p_true|). Per pT bin and parameter the figure
// shows the mean of the residuals with their standard deviation as error bar, and as
// a band around zero +- the median sqrt(cov) of the estimate (err_loc0, err_loc1,
// err_phi, err_theta, err_qop / |q/p|). A mean away from zero is a bias; an error bar
// wider than the band means the covariance is too confident.
//
// Output in figures/seed_parameters/<input_dir>/, with
// <name> = seed_parameter_resolution_iter_<iter>_<own_hits>of3{,_seedexit,_secondpassstart}
// (no suffix for pass 0):
//  - <name>.pdf, and <name>.root with the per-bin mean, standard deviation and median
//    sqrt(cov) (h_mean_*, h_sigma_*),
//  - <name>_distributions.root: the residual distribution of every parameter and pT
//    bin (directory per parameter, h_<parameter>_bin<i>), and the tree "residuals"
//    with one entry per used seed (event, seed index in the event, particle barcode,
//    true pT, residuals and sqrt(cov)), to look up single particles.
//
// Usage (from low_pt_studies):
//   root -l -b -q 'macros/seed_parameter_resolution.C+("gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_OT_seeding_matching_2clus", 4, 3)'
//   root -l -b -q 'macros/seed_parameter_resolution.C+("gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_OT_seeding_matching_2clus_initial_sigmas", 4, 3, 1)'
namespace {

double median(std::vector<double> values) {
  if (values.empty()) return 0;
  const size_t index = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + index, values.end());
  return values[index];
}

}  // namespace

void seed_parameter_resolution(const TString input_dir = "gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_OT_seeding_matching_2clus",
                               const int iteration = 4,
                               const int own_hits = 3,
                               const int pass = 0,
                               const bool onSTBC = false) {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);

  if (pass < 0 || pass > 2) {
    std::cout << "Error: pass must be 0 (seed estimate), 1 (after the seed) or 2 (second pass start)"
              << std::endl;
    return;
  }
  const bool exit_states = pass > 0;
  const TString prefix = iteration > 0 ? Form("iter_%d_", iteration) : "";
  const TString filename = Utils::outputPathBase(onSTBC) + input_dir + Form("/iter_%d/", iteration) + prefix +
                           (exit_states ? "seedexitparameters.root" : "estimatedparameters.root");
  const char* tree_name = exit_states ? "seedexitparams" : "estimatedparams";
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << filename << std::endl;
    return;
  }
  TTree* tree = (TTree*)file->Get(tree_name);
  if (!tree) {
    std::cout << "Error: " << tree_name << " tree not found in " << filename << std::endl;
    return;
  }

  Bool_t truth_matched;
  Float_t t_pt, qop, t_qop, phi, t_phi, theta, t_theta, err_qop, err_phi, err_theta, res_loc0, res_loc1,
      err_loc0, err_loc1;
  tree->SetBranchAddress("truthMatched", &truth_matched);
  tree->SetBranchAddress("t_pt", &t_pt);
  tree->SetBranchAddress("qop", &qop);
  tree->SetBranchAddress("t_qop", &t_qop);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("t_phi", &t_phi);
  tree->SetBranchAddress("theta", &theta);
  tree->SetBranchAddress("t_theta", &t_theta);
  tree->SetBranchAddress("err_qop", &err_qop);
  tree->SetBranchAddress("err_phi", &err_phi);
  tree->SetBranchAddress("err_theta", &err_theta);
  tree->SetBranchAddress("res_loc0", &res_loc0);
  tree->SetBranchAddress("res_loc1", &res_loc1);
  tree->SetBranchAddress("err_loc0", &err_loc0);
  tree->SetBranchAddress("err_loc1", &err_loc1);
  // seed exit states: entry -> seed, pass and whether the state exists
  UInt_t exit_event = 0;
  Int_t exit_seed_index = 0, exit_pass = 0, exit_status = 0;
  Bool_t t_on_particle = true;
  Int_t seed_space_point = -1, n_seed_hits_on_track = -1;
  if (exit_states) {
    tree->SetBranchAddress("seed_space_point", &seed_space_point);
    tree->SetBranchAddress("n_seed_hits_on_track", &n_seed_hits_on_track);
    tree->SetBranchAddress("event_nr", &exit_event);
    tree->SetBranchAddress("seed_index", &exit_seed_index);
    tree->SetBranchAddress("pass", &exit_pass);
    tree->SetBranchAddress("status", &exit_status);
    tree->SetBranchAddress("t_on_particle", &t_on_particle);
  }

  // the triplet of every seed, in the order of the parameters
  const TString iter_dir = Utils::outputPathBase(onSTBC) + input_dir + Form("/iter_%d/", iteration);
  TFile* seed_file = TFile::Open(iter_dir + "seeds.root");
  TTree* seed_tree = seed_file ? (TTree*)seed_file->Get("seeds") : nullptr;
  // the estimated parameters are aligned entry by entry, the exit states (two per
  // seed) by event and seed index
  if (!seed_tree || (!exit_states && seed_tree->GetEntries() != tree->GetEntries()) ||
      (exit_states && 2 * seed_tree->GetEntries() != tree->GetEntries())) {
    std::cout << "Error: seeds.root missing or not aligned with the parameters" << std::endl;
    return;
  }
  UInt_t seed_event;
  ULong64_t seed_measurement[3];
  seed_tree->SetBranchAddress("event_id", &seed_event);
  for (int i_sp = 0; i_sp < 3; i_sp++) {
    seed_tree->SetBranchAddress(Form("measurement_id_%d", i_sp + 1), &seed_measurement[i_sp]);
  }
  using Key = std::pair<long, unsigned long>;  // event, measurement index
  std::vector<std::array<Key, 3>> triplets(seed_tree->GetEntries());
  std::set<Key> needed;
  std::map<std::pair<long, int>, Long64_t> seed_entry;  // (event, seed index in event) -> entry
  std::map<long, int> n_event_seeds;
  std::vector<int> seed_index_in_event(seed_tree->GetEntries());
  for (Long64_t i_entry = 0; i_entry < seed_tree->GetEntries(); i_entry++) {
    seed_tree->GetEntry(i_entry);
    seed_index_in_event[i_entry] = n_event_seeds[seed_event];
    seed_entry[{seed_event, n_event_seeds[seed_event]++}] = i_entry;
    for (int i_sp = 0; i_sp < 3; i_sp++) {
      triplets[i_entry][i_sp] = {seed_event, seed_measurement[i_sp]};
      needed.insert(triplets[i_entry][i_sp]);
    }
  }

  // truth particles of those measurements
  using Barcode = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned>;
  std::map<Key, std::vector<Barcode>> measurement_particles;
  {
    TFile* measurement_file = TFile::Open(Utils::outputPathBase(onSTBC) + input_dir + "/measurements_with_noise.root");
    TTree* measurement_tree = measurement_file ? (TTree*)measurement_file->Get("measurements") : nullptr;
    if (!measurement_tree) {
      std::cout << "Error: measurements_with_noise.root missing" << std::endl;
      return;
    }
    Int_t measurement_event;
    ULong_t measurement_id;
    std::vector<unsigned>*vertex_primary = nullptr, *vertex_secondary = nullptr, *particle = nullptr,
    *generation = nullptr, *sub_particle = nullptr;
    measurement_tree->SetBranchStatus("*", false);
    for (const char* name : {"event_nr", "measurement_id", "particles_vertex_primary", "particles_vertex_secondary",
                             "particles_particle", "particles_generation", "particles_sub_particle"}) {
      measurement_tree->SetBranchStatus(name, true);
    }
    measurement_tree->SetBranchAddress("event_nr", &measurement_event);
    measurement_tree->SetBranchAddress("measurement_id", &measurement_id);
    measurement_tree->SetBranchAddress("particles_vertex_primary", &vertex_primary);
    measurement_tree->SetBranchAddress("particles_vertex_secondary", &vertex_secondary);
    measurement_tree->SetBranchAddress("particles_particle", &particle);
    measurement_tree->SetBranchAddress("particles_generation", &generation);
    measurement_tree->SetBranchAddress("particles_sub_particle", &sub_particle);
    // only the events with seeds
    std::set<long> seed_events;
    for (const auto& key : needed) seed_events.insert(key.first);
    TBranch* event_branch = measurement_tree->GetBranch("event_nr");
    for (Long64_t i_entry = 0; i_entry < measurement_tree->GetEntries(); i_entry++) {
      event_branch->GetEntry(i_entry);
      if (!seed_events.count(measurement_event)) continue;
      measurement_tree->GetEntry(i_entry);
      const Key key{measurement_event, measurement_id};
      if (!needed.count(key)) continue;
      auto& barcodes = measurement_particles[key];
      for (size_t i_particle = 0; i_particle < particle->size(); i_particle++) {
        barcodes.emplace_back((*vertex_primary)[i_particle], (*vertex_secondary)[i_particle],
                              (*particle)[i_particle], (*generation)[i_particle], (*sub_particle)[i_particle]);
      }
    }
    measurement_file->Close();
  }

  UInt_t majority_vertex_primary, majority_vertex_secondary, majority_particle, majority_generation, majority_sub_particle;
  tree->SetBranchAddress("particleId_vertex_primary", &majority_vertex_primary);
  tree->SetBranchAddress("particleId_vertex_secondary", &majority_vertex_secondary);
  tree->SetBranchAddress("particleId_particle", &majority_particle);
  tree->SetBranchAddress("particleId_generation", &majority_generation);
  tree->SetBranchAddress("particleId_sub_particle", &majority_sub_particle);

  // true pT bins in MeV: <25, 25-50, 50-75, 75-100, 100-150, 150-200, >200 (drawn up to 500)
  const std::vector<double> edges = {0, 25, 50, 75, 100, 150, 200, 500};
  const int n_bins = edges.size() - 1;
  const double same_surface_tolerance = 0.1;  // mm

  // [parameter][bin]: residuals and sqrt(cov) of the used seeds;
  // parameters loc0, loc1 (um), phi, theta (deg), q/p (relative)
  const int n_parameters = 5;
  std::vector<std::vector<double>> residuals[n_parameters], sigmas[n_parameters];
  for (int i_par = 0; i_par < n_parameters; i_par++) {
    residuals[i_par].resize(n_bins);
    sigmas[i_par].resize(n_bins);
  }

  // one entry per used seed, for the tree of the distributions file
  struct SeedRecord {
    long event;
    int seed_index;
    unsigned barcode[5];
    double t_pt;
    int pt_bin;
    int seed_space_point;
    double residual[n_parameters];
    double sigma[n_parameters];
  };
  std::vector<SeedRecord> records;

  long n_seeds = 0, n_unmatched = 0, n_few_own_hits = 0, n_other_surface = 0, n_opposite_charge = 0, n_clean = 0;
  long n_no_state = 0;  // exit states: no state of this pass (status != 0)
  long n_short = 0;     // pass 1: the CKF did not take the seed up to its end
  for (Long64_t i_entry = 0; i_entry < tree->GetEntries(); i_entry++) {
    tree->GetEntry(i_entry);
    Long64_t i_seed = i_entry;
    if (exit_states) {
      if (exit_pass != pass) continue;
      const auto it = seed_entry.find({long(exit_event), exit_seed_index});
      if (it == seed_entry.end()) {
        std::cout << "Error: seed " << exit_seed_index << " of event " << exit_event << " not in seeds.root"
                  << std::endl;
        return;
      }
      i_seed = it->second;
    }
    n_seeds++;
    if (!truth_matched) {
      n_unmatched++;
      continue;
    }
    if (exit_states && exit_status != 0) {
      n_no_state++;
      continue;
    }
    // pass 1: the state has to be at the end of the seed, i.e. on its third space
    // point or its TOF hit, with the three tracker hits all taken by the CKF
    if (pass == 1 && (seed_space_point < 2 || n_seed_hits_on_track < 3)) {
      n_short++;
      continue;
    }
    const Barcode majority{majority_vertex_primary, majority_vertex_secondary, majority_particle,
                           majority_generation, majority_sub_particle};
    int n_own = 0;
    for (const auto& key : triplets[i_seed]) {
      const auto it = measurement_particles.find(key);
      if (it != measurement_particles.end() &&
          std::find(it->second.begin(), it->second.end(), majority) != it->second.end()) {
        n_own++;
      }
    }
    if (n_own != own_hits) {
      n_few_own_hits++;
      continue;
    }
    // pass 0: truth and parameters on different surfaces; pass 1, 2: the state sits on
    // a hit without a hit of the particle
    if (exit_states ? !t_on_particle
                    : (std::abs(res_loc0) > same_surface_tolerance || std::abs(res_loc1) > same_surface_tolerance)) {
      n_other_surface++;
      continue;
    }
    if (qop * t_qop < 0) {
      n_opposite_charge++;
      continue;
    }
    n_clean++;
    const double pt_mev = 1000. * t_pt;  // ACTS units: GeV
    const int i_bin = std::upper_bound(edges.begin() + 1, edges.end() - 1, pt_mev) - (edges.begin() + 1);
    const double mm_to_um = 1000.;
    residuals[0][i_bin].push_back(mm_to_um * res_loc0);
    sigmas[0][i_bin].push_back(mm_to_um * err_loc0);
    residuals[1][i_bin].push_back(mm_to_um * res_loc1);
    sigmas[1][i_bin].push_back(mm_to_um * err_loc1);
    residuals[2][i_bin].push_back(TMath::RadToDeg() * TVector2::Phi_mpi_pi(phi - t_phi));
    sigmas[2][i_bin].push_back(TMath::RadToDeg() * err_phi);
    residuals[3][i_bin].push_back(TMath::RadToDeg() * (theta - t_theta));
    sigmas[3][i_bin].push_back(TMath::RadToDeg() * err_theta);
    residuals[4][i_bin].push_back((std::abs(qop) - std::abs(t_qop)) / std::abs(t_qop));
    sigmas[4][i_bin].push_back(err_qop / std::abs(qop));

    SeedRecord record{triplets[i_seed][0].first,
                      exit_states ? exit_seed_index : seed_index_in_event[i_seed],
                      {majority_vertex_primary, majority_vertex_secondary, majority_particle,
                       majority_generation, majority_sub_particle},
                      pt_mev,
                      i_bin,
                      exit_states ? seed_space_point : -1,
                      {},
                      {}};
    for (int i_par = 0; i_par < n_parameters; i_par++) {
      record.residual[i_par] = residuals[i_par][i_bin].back();
      record.sigma[i_par] = sigmas[i_par][i_bin].back();
    }
    records.push_back(record);
  }
  std::cout << "Seeds: " << n_seeds << ", used " << n_clean << " (not truth-matched " << n_unmatched
            << ", not " << own_hits << "/3 own hits " << n_few_own_hits << ", truth on another surface "
            << n_other_surface << ", opposite charge " << n_opposite_charge;
  if (exit_states) std::cout << ", no state of pass " << pass << " " << n_no_state;
  if (pass == 1) std::cout << ", seed not taken to its end " << n_short;
  std::cout << ")" << std::endl;

  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/low_pt_studies/figures/seed_parameters/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);
  const char* pass_suffix[3] = {"", "_seedexit", "_secondpassstart"};
  const TString output_name =
      output_dir + Form("/seed_parameter_resolution_iter_%d_%dof3%s", iteration, own_hits, pass_suffix[pass]);
  TFile* output = TFile::Open(output_name + ".root", "RECREATE");

  const char* names[n_parameters] = {"loc0", "loc1", "phi", "theta", "qop"};
  const char* titles[n_parameters] = {"loc_{0} [#mum]", "loc_{1} [#mum]", "#phi [deg]", "#theta [deg]",
                                      "(|q/p| - |q/p_{true}|) / |q/p_{true}|"};
  TCanvas* canvas = new TCanvas("c_seed_parameter_resolution",
                                Form("Seed parameter resolution, iteration %d", iteration), 2400, 1200);
  canvas->Divide(3, 2);
  for (int i_par = 0; i_par < n_parameters; i_par++) {
    // mean of the signed residuals, standard deviation as error bar
    TH1D* mean = new TH1D(Form("h_mean_%s", names[i_par]),
                          Form(";true p_{T} [MeV];%s, estimate - truth", titles[i_par]), n_bins, edges.data());
    // median sqrt(cov), drawn as +- band around zero
    TH1D* sigma = new TH1D(Form("h_sigma_%s", names[i_par]), Form(";true p_{T} [MeV];%s", titles[i_par]),
                           n_bins, edges.data());
    double range = 0;
    for (int i_bin = 0; i_bin < n_bins; i_bin++) {
      const auto& values = residuals[i_par][i_bin];
      if (values.empty()) continue;
      double sum = 0, sum2 = 0;
      for (double value : values) sum += value;
      const double average = sum / values.size();
      for (double value : values) sum2 += (value - average) * (value - average);
      const double deviation = values.size() > 1 ? std::sqrt(sum2 / (values.size() - 1)) : 0;
      mean->SetBinContent(i_bin + 1, average);
      mean->SetBinError(i_bin + 1, deviation);
      sigma->SetBinContent(i_bin + 1, median(sigmas[i_par][i_bin]));
      range = std::max({range, std::abs(average) + deviation, sigma->GetBinContent(i_bin + 1)});
    }
    output->cd();
    mean->Write();
    sigma->Write();

    TH1D* band = (TH1D*)sigma->Clone(Form("h_band_%s", names[i_par]));
    for (int i_bin = 1; i_bin <= n_bins; i_bin++) {
      band->SetBinError(i_bin, band->GetBinContent(i_bin));
      band->SetBinContent(i_bin, 0);
    }
    band->SetFillColorAlpha(kRed + 1, 0.25);
    band->SetLineColor(kRed + 1);
    band->SetMarkerSize(0);
    mean->SetMarkerStyle(20);
    mean->SetMarkerColor(kBlue + 1);
    mean->SetLineColor(kBlue + 1);
    mean->SetLineWidth(2);

    canvas->cd(i_par + 1);
    // symmetric around zero, 10% beyond the largest |mean +- std. dev.| or sqrt(cov)
    mean->SetMaximum(1.1 * range);
    mean->SetMinimum(-1.1 * range);
    mean->Draw("E1");
    band->Draw("E2 same");
    mean->Draw("E1 same");
    TLine* zero = new TLine(edges.front(), 0, edges.back(), 0);
    zero->SetLineStyle(2);
    zero->Draw();
    mean->GetXaxis()->SetNdivisions(510);

    TLegend* legend = new TLegend(0.40, 0.78, 0.88, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(mean, exit_states ? "CKF - truth: mean #pm std. dev." : "seed estimate - truth: mean #pm std. dev.",
                     "lep");
    legend->AddEntry(band, exit_states ? "#pm CKF #sqrt{cov} (median)" : "#pm initial #sigma (median)", "f");
    legend->Draw();
  }

  canvas->cd(6);
  TPaveText* summary = new TPaveText(0.05, 0.1, 0.95, 0.9, "NDC");
  summary->SetBorderSize(0);
  summary->SetFillStyle(0);
  summary->SetTextAlign(12);
  const char* pass_title[3] = {"seed estimate (CKF start)", "after the seed (filtered on its last hit)",
                               "second pass start (smoothed on the first hit)"};
  summary->AddText(Form("iteration %d, %s: %ld seeds", iteration, pass_title[pass], n_seeds));
  summary->AddText(Form("used (%d/3 hits of the particle): %ld", own_hits, n_clean));
  summary->AddText(Form("not truth-matched: %ld", n_unmatched));
  summary->AddText(Form("other number of hits of the particle: %ld", n_few_own_hits));
  if (exit_states) {
    summary->AddText(Form("no state of this pass: %ld", n_no_state));
    if (pass == 1) summary->AddText(Form("CKF did not take the seed to its end: %ld", n_short));
    summary->AddText(Form("state on a hit of another particle: %ld", n_other_surface));
  } else {
    summary->AddText(Form("truth on another surface: %ld", n_other_surface));
  }
  summary->AddText(Form("opposite charge to truth: %ld", n_opposite_charge));
  TString counts = "clean seeds per bin:";
  for (int i_bin = 0; i_bin < n_bins; i_bin++) counts += Form(" %zu", residuals[0][i_bin].size());
  summary->AddText(counts);
  summary->AddText("last bin: p_{T} > 200 MeV");
  summary->Draw();

  canvas->SaveAs(output_name + ".pdf");
  output->Close();
  std::cout << "Wrote " << output_name << ".{pdf,root}" << std::endl;

  // the residual distributions per parameter and pT bin, and one entry per used seed
  TFile* distributions = TFile::Open(output_name + "_distributions.root", "RECREATE");
  const int n_histogram_bins = 40;
  for (int i_par = 0; i_par < n_parameters; i_par++) {
    distributions->mkdir(names[i_par])->cd();
    for (int i_bin = 0; i_bin < n_bins; i_bin++) {
      const auto& values = residuals[i_par][i_bin];
      double low = -1, high = 1;
      if (!values.empty()) {
        const auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
        const double pad = std::max(0.05 * (*max_it - *min_it), 1e-6 * std::max(1., std::abs(*max_it)));
        low = *min_it - pad;
        high = *max_it + pad;
      }
      const TString pt_range = i_bin == 0           ? Form("< %.0f MeV", edges[1])
                               : i_bin + 1 < n_bins ? Form("%.0f-%.0f MeV", edges[i_bin], edges[i_bin + 1])
                                                    : Form("> %.0f MeV", edges[i_bin]);
      TH1D* histogram = new TH1D(Form("h_%s_bin%d", names[i_par], i_bin),
                                 Form("true p_{T} %s;%s, estimate - truth;seeds", pt_range.Data(), titles[i_par]),
                                 n_histogram_bins, low, high);
      for (double value : values) histogram->Fill(value);
      histogram->Write();
    }
  }
  distributions->cd();
  TTree* residual_tree = new TTree("residuals", "residuals of the used seeds");
  Long64_t out_event;
  Int_t out_seed_index, out_pt_bin, out_seed_space_point;
  UInt_t out_barcode[5];
  Double_t out_t_pt, out_residual[n_parameters], out_sigma[n_parameters];
  residual_tree->Branch("event", &out_event);
  residual_tree->Branch("seed_index", &out_seed_index);
  const char* barcode_names[5] = {"vertex_primary", "vertex_secondary", "particle", "generation", "sub_particle"};
  for (int i_id = 0; i_id < 5; i_id++) residual_tree->Branch(barcode_names[i_id], &out_barcode[i_id]);
  residual_tree->Branch("t_pt", &out_t_pt);  // MeV
  residual_tree->Branch("pt_bin", &out_pt_bin);
  residual_tree->Branch("seed_space_point", &out_seed_space_point);  // pass 1, 2: hit of the state
  for (int i_par = 0; i_par < n_parameters; i_par++) {
    residual_tree->Branch(Form("res_%s", names[i_par]), &out_residual[i_par]);
    residual_tree->Branch(Form("sigma_%s", names[i_par]), &out_sigma[i_par]);
  }
  for (const SeedRecord& record : records) {
    out_event = record.event;
    out_seed_index = record.seed_index;
    for (int i_id = 0; i_id < 5; i_id++) out_barcode[i_id] = record.barcode[i_id];
    out_t_pt = record.t_pt;
    out_pt_bin = record.pt_bin;
    out_seed_space_point = record.seed_space_point;
    for (int i_par = 0; i_par < n_parameters; i_par++) {
      out_residual[i_par] = record.residual[i_par];
      out_sigma[i_par] = record.sigma[i_par];
    }
    residual_tree->Fill();
  }
  residual_tree->Write();
  distributions->Close();
  std::cout << "Wrote " << output_name << "_distributions.root" << std::endl;
}
