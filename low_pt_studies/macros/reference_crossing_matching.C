#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../../Utils/BasicFunctions.h"

// Differences between tracks at the reference cylinder (reference_crossings.root,
// RootReferenceCrossingWriter in ACTSO2), as the basis for matching the track segments
// of one particle (e.g. the legs of a looper).
//
// Every pair of tracks in an event is compared once: a looper crosses the cylinder
// several times, so of all combinations of a crossing of one track with a crossing of
// the other only the closest one in 3D is used (purely spatial). The pairs are split by
// truth: both tracks with the same majority particle ("same particle"), and all others.
// Each panel overlays both, normalised to unit area. The time difference is shown
// separately for pairs where both tracks have a TOF hit (iTOF or oTOF) and for the
// others, as is its pull (time difference over the combined time uncertainty).
//
// Histogram ranges: from the 0.5% and 99.5% quantiles of both classes together
// (symmetric for signed differences), so no range is tuned by hand.
//
// Output: figures/reference_crossings/<input_dir>/reference_crossing_matching.{pdf,root};
// the ROOT file holds the canvas and all histograms.
//
// Usage (from low_pt_studies):
//   root -l -b -q 'macros/reference_crossing_matching.C("gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_2")'
namespace {

struct Crossing {
  float x, y, z, phi, dirPhi, theta, qop, time, timeSigma, helixRadius;
  int radialSign;
};

struct Track {
  unsigned trackNr;
  unsigned long long particleHash;
  bool hasTof;
  std::vector<Crossing> crossings;
};

// angle difference wrapped to [-pi, pi)
double deltaAngle(double a, double b) { return std::remainder(a - b, 2 * M_PI); }

// one quantity, filled for both truth classes
struct Quantity {
  TString name, title;
  bool symmetric;
  std::vector<double> values[2] = {};  // 0: same particle, 1: other
};

double quantile(std::vector<double> values, double q) {
  if (values.empty()) return 0;
  const size_t index = std::min(values.size() - 1, size_t(q * (values.size() - 1)));
  std::nth_element(values.begin(), values.begin() + index, values.end());
  return values[index];
}

}  // namespace

void reference_crossing_matching(const TString input_dir = "gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_2",
                                 const bool onSTBC = false,
                                 const int n_bins = 100) {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);

  const TString filename = Utils::outputPathBase(onSTBC) + input_dir + "/reference_crossings.root";
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << filename << std::endl;
    return;
  }
  TTree* tree = (TTree*)file->Get("crossings");
  if (!tree) {
    std::cout << "Error: crossings tree not found in " << filename << std::endl;
    return;
  }

  unsigned event_nr, track_nr, n_itof_hits, n_otof_hits;
  float radius;
  ULong_t t_particle_hash = 0;  // std::uint64_t in the writer
  std::vector<float>*x = nullptr, *y = nullptr, *z = nullptr, *phi = nullptr, *dir_phi = nullptr,
  *theta = nullptr, *qop = nullptr, *time = nullptr, *time_sigma = nullptr, *helix_radius = nullptr;
  std::vector<int>* radial_sign = nullptr;
  tree->SetBranchAddress("event_nr", &event_nr);
  tree->SetBranchAddress("track_nr", &track_nr);
  tree->SetBranchAddress("radius", &radius);
  tree->SetBranchAddress("n_itof_hits", &n_itof_hits);
  tree->SetBranchAddress("n_otof_hits", &n_otof_hits);
  tree->SetBranchAddress("t_particle_hash", &t_particle_hash);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("dir_phi", &dir_phi);
  tree->SetBranchAddress("theta", &theta);
  tree->SetBranchAddress("qop", &qop);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("time_sigma", &time_sigma);
  tree->SetBranchAddress("helix_radius", &helix_radius);
  tree->SetBranchAddress("radial_sign", &radial_sign);

  // tracks with at least one crossing, per event (entries of an event need not be
  // contiguous when the chain ran multithreaded)
  std::map<unsigned, std::vector<Track>> events;
  float reference_radius = 0;
  for (Long64_t i_entry = 0; i_entry < tree->GetEntries(); i_entry++) {
    tree->GetEntry(i_entry);
    reference_radius = radius;
    if (x->empty()) continue;
    Track track{track_nr, t_particle_hash, n_itof_hits + n_otof_hits > 0, {}};
    for (size_t i_crossing = 0; i_crossing < x->size(); i_crossing++) {
      track.crossings.push_back({(*x)[i_crossing], (*y)[i_crossing], (*z)[i_crossing],
                                 (*phi)[i_crossing], (*dir_phi)[i_crossing], (*theta)[i_crossing],
                                 (*qop)[i_crossing], (*time)[i_crossing], (*time_sigma)[i_crossing],
                                 (*helix_radius)[i_crossing], (*radial_sign)[i_crossing]});
    }
    events[event_nr].push_back(std::move(track));
  }

  std::vector<Quantity> quantities = {
      {"distance", "closest crossings: |#Deltaxyz| [mm]", false},
      {"rdphi", "R #Delta#phi [mm]", true},
      {"dz", "#Deltaz [mm]", true},
      {"ddirphi", "#Delta#phi_{direction} [rad]", true},
      {"dtheta", "#Delta#theta [rad]", true},
      {"radial", "radial direction: product of the signs (+1 same, -1 opposite)", true},
      {"charge", "charge: product of the signs of q/p (+1 same, -1 opposite)", true},
      {"drho", "helix radius: 2 (#rho_{1} - #rho_{2}) / (#rho_{1} + #rho_{2})", true},
      {"dt_tof", "#Deltat, both tracks with a TOF hit [mm/c]", true},
      {"dt_notof", "#Deltat, not both with a TOF hit [mm/c]", true},
      {"pull_tof", "#Deltat / #sigma_{#Deltat}, both with a TOF hit", true},
      {"pull_notof", "#Deltat / #sigma_{#Deltat}, not both with a TOF hit", true},
  };
  auto fill = [&](const TString& name, int truth_class, double value) {
    for (auto& quantity : quantities) {
      if (quantity.name == name) {
        if (std::isfinite(value)) quantity.values[truth_class].push_back(value);
        return;
      }
    }
  };

  size_t n_pairs[2] = {0, 0};
  for (const auto& [event, tracks] : events) {
    for (size_t i_track = 0; i_track < tracks.size(); i_track++) {
      for (size_t j_track = i_track + 1; j_track < tracks.size(); j_track++) {
        const Track& first = tracks[i_track];
        const Track& second = tracks[j_track];
        // closest pair of crossings in 3D
        double best_distance = 1e30;
        const Crossing* a = nullptr;
        const Crossing* b = nullptr;
        for (const auto& crossing_a : first.crossings) {
          for (const auto& crossing_b : second.crossings) {
            const double distance = std::sqrt(std::pow(crossing_a.x - crossing_b.x, 2) +
                                              std::pow(crossing_a.y - crossing_b.y, 2) +
                                              std::pow(crossing_a.z - crossing_b.z, 2));
            if (distance < best_distance) {
              best_distance = distance;
              a = &crossing_a;
              b = &crossing_b;
            }
          }
        }
        const int truth_class =
            (first.particleHash != 0 && first.particleHash == second.particleHash) ? 0 : 1;
        n_pairs[truth_class]++;

        fill("distance", truth_class, best_distance);
        fill("rdphi", truth_class, reference_radius * deltaAngle(a->phi, b->phi));
        fill("dz", truth_class, a->z - b->z);
        fill("ddirphi", truth_class, deltaAngle(a->dirPhi, b->dirPhi));
        fill("dtheta", truth_class, a->theta - b->theta);
        fill("radial", truth_class, a->radialSign * b->radialSign);
        fill("charge", truth_class, (a->qop * b->qop > 0) ? 1 : -1);
        fill("drho", truth_class, 2 * (a->helixRadius - b->helixRadius) / (a->helixRadius + b->helixRadius));
        const bool both_tof = first.hasTof && second.hasTof;
        const double dt = a->time - b->time;
        const double sigma_dt = std::hypot(a->timeSigma, b->timeSigma);
        fill(both_tof ? "dt_tof" : "dt_notof", truth_class, dt);
        fill(both_tof ? "pull_tof" : "pull_notof", truth_class, sigma_dt > 0 ? dt / sigma_dt : NAN);
      }
    }
  }
  std::cout << "Reference cylinder r = " << reference_radius << " mm, " << events.size()
            << " events with crossings; track pairs: " << n_pairs[0] << " same particle, " << n_pairs[1]
            << " other" << std::endl;

  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/low_pt_studies/figures/reference_crossings/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);
  TFile* output = TFile::Open(output_dir + "/reference_crossing_matching.root", "RECREATE");

  TCanvas* canvas = new TCanvas("c_reference_crossing_matching",
                                Form("Track pairs at the reference cylinder r = %.0f mm", reference_radius),
                                2400, 1500);
  canvas->Divide(4, 3);
  const char* class_names[2] = {"same particle", "other"};
  const int class_colors[2] = {kBlue + 1, kRed + 1};
  for (size_t i_quantity = 0; i_quantity < quantities.size(); i_quantity++) {
    Quantity& quantity = quantities[i_quantity];
    std::vector<double> all = quantity.values[0];
    all.insert(all.end(), quantity.values[1].begin(), quantity.values[1].end());
    double low = quantile(all, 0.005), high = quantile(all, 0.995);
    if (quantity.symmetric) {
      high = std::max(std::abs(low), std::abs(high));
      low = -high;
    }
    // discrete sign products: one bin per value
    const bool is_sign = quantity.name == "radial" || quantity.name == "charge";
    const int bins = is_sign ? 3 : n_bins;
    if (is_sign) {
      low = -1.5;
      high = 1.5;
    }
    if (!(high > low)) high = low + 1;

    canvas->cd(i_quantity + 1);
    gPad->SetLogy();
    TLegend* legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    double maximum = 0;
    TH1D* histograms[2];
    for (int truth_class = 0; truth_class < 2; truth_class++) {
      TH1D* histogram = new TH1D(Form("h_%s_%s", quantity.name.Data(), truth_class == 0 ? "same" : "other"),
                                 Form(";%s;fraction of pairs", quantity.title.Data()), bins, low, high);
      for (double value : quantity.values[truth_class]) histogram->Fill(value);
      histogram->SetLineColor(class_colors[truth_class]);
      histogram->SetLineWidth(2);
      output->cd();
      histogram->Write();
      if (histogram->Integral(0, bins + 1) > 0) histogram->Scale(1. / histogram->Integral(0, bins + 1));
      maximum = std::max(maximum, histogram->GetMaximum());
      legend->AddEntry(histogram, Form("%s (%zu)", class_names[truth_class], quantity.values[truth_class].size()), "l");
      histograms[truth_class] = histogram;
    }
    for (int truth_class = 0; truth_class < 2; truth_class++) {
      histograms[truth_class]->SetMaximum(maximum * 3);
      histograms[truth_class]->SetMinimum(maximum * 1e-4);
      histograms[truth_class]->Draw(truth_class == 0 ? "hist" : "hist same");
    }
    legend->Draw();
  }
  canvas->SaveAs(output_dir + "/reference_crossing_matching.pdf");
  output->cd();
  canvas->Write();
  output->Close();
  std::cout << "Wrote " << output_dir << "/reference_crossing_matching.{pdf,root}" << std::endl;
}
