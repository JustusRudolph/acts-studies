#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <Math/PdfFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../../Utils/BasicFunctions.h"

// chi2 of the track joins (track_joins.root, RootTrackJoinWriter in ACTSO2): the
// earlier track a is propagated from its last measurement to the curvilinear plane at
// the first measurement of the later track b, and both are compared there with their
// summed covariances (loc0, loc1, phi, theta, q/p, t).
//
// Pairs are split by truth: both tracks with the same (nonzero) majority particle
// ("same particle"), and all others. Only joins whose propagation reached the plane
// (ok) enter the distributions; the counts of failed ones are in the last panel.
//
// Panels: log10 chi2 without the time (5 dof), with the expected chi2(5) density for
// reference; log10 chi2 with the time (6 dof) for pairs where both tracks have a TOF
// hit and for the others; the p-value of chi2 (flat for same-particle pairs if the
// covariances are right); the pull of each parameter, delta / sqrt(sigma_a^2 +
// sigma_b^2) (diagonal only, unit width if right); log10 chi2 of the self joins (a
// track's last measurement propagated back to its own first one, if written) against
// the same-particle pairs. Self joins are not part of the two truth classes. The
// same-particle chi2 per direction of the propagation is in the summary panel.
//
// Histogram ranges: from the 0.5% and 99.5% quantiles of both classes together
// (symmetric for pulls), so no range is tuned by hand.
//
// Output: figures/track_joins/<input_dir>/track_join_chi2.{pdf,root}.
//
// Usage (from low_pt_studies):
//   root -l -b -q 'macros/track_join_chi2.C("gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_OT_seeding_matching_1")'
namespace {

// one quantity, filled for two classes
struct Quantity {
  TString name, title;
  bool symmetric;
  TString labels[2] = {"same particle", "other"};
  std::vector<double> values[2] = {};
};

double quantile(std::vector<double> values, double q) {
  if (values.empty()) return 0;
  const size_t index = std::min(values.size() - 1, size_t(q * (values.size() - 1)));
  std::nth_element(values.begin(), values.begin() + index, values.end());
  return values[index];
}

}  // namespace

void track_join_chi2(const TString input_dir = "gamma_conversion/seed_1_reco_ortho_seed_chord_navigator_OT_seeding_matching_1",
                     const bool onSTBC = false,
                     const int n_bins = 60) {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);

  const TString filename = Utils::outputPathBase(onSTBC) + input_dir + "/track_joins.root";
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << filename << std::endl;
    return;
  }
  TTree* tree = (TTree*)file->Get("joins");
  if (!tree) {
    std::cout << "Error: joins tree not found in " << filename << std::endl;
    return;
  }

  ULong_t hash_a = 0, hash_b = 0;  // std::uint64_t in the writer
  unsigned n_itof_a, n_otof_a, n_itof_b, n_otof_b;
  float sigma_a[6], sigma_b[6], delta[6], chi2, chi2_time;
  bool ok, backward, self = false;
  tree->SetBranchAddress("t_particle_hash_a", &hash_a);
  tree->SetBranchAddress("t_particle_hash_b", &hash_b);
  tree->SetBranchAddress("n_itof_hits_a", &n_itof_a);
  tree->SetBranchAddress("n_otof_hits_a", &n_otof_a);
  tree->SetBranchAddress("n_itof_hits_b", &n_itof_b);
  tree->SetBranchAddress("n_otof_hits_b", &n_otof_b);
  tree->SetBranchAddress("sigma_a", sigma_a);
  tree->SetBranchAddress("sigma_b", sigma_b);
  tree->SetBranchAddress("delta", delta);
  tree->SetBranchAddress("chi2", &chi2);
  tree->SetBranchAddress("chi2_time", &chi2_time);
  tree->SetBranchAddress("ok", &ok);
  tree->SetBranchAddress("backward", &backward);
  const bool has_self = tree->GetBranch("self") != nullptr;
  if (has_self) tree->SetBranchAddress("self", &self);
  // written since TrackJoiningAlgorithm: chi2 with q/p free (4 dof) and |p_b| / |p_a|
  float chi2_free_qop = -1, momentum_ratio = -1;
  const bool has_free_qop = tree->GetBranch("chi2_free_qop") != nullptr;
  if (has_free_qop) {
    tree->SetBranchAddress("chi2_free_qop", &chi2_free_qop);
    tree->SetBranchAddress("momentum_ratio", &momentum_ratio);
  }

  const char* parameter_names[6] = {"loc0", "loc1", "#phi", "#theta", "q/p", "t"};
  std::vector<Quantity> quantities = {
      {"chi2", "log_{10} #chi^{2} (5 dof, no time)", false},
      {"chi2_time_tof", "log_{10} #chi^{2} (6 dof), both tracks with TOF hit", false},
      {"chi2_time_notof", "log_{10} #chi^{2} (6 dof), not both with TOF hit", false},
      {"prob", "p-value of #chi^{2} (5 dof)", false},
  };
  for (int i_par = 0; i_par < 6; i_par++) {
    quantities.push_back({Form("pull_%d", i_par), Form("pull of %s", parameter_names[i_par]), true});
  }
  quantities.push_back({"chi2_self", "log_{10} #chi^{2} (5 dof)", false, {"self join", "same particle"}});
  if (has_free_qop) {
    quantities.push_back({"chi2_free_qop", "log_{10} #chi^{2} with q/p free (4 dof)", false});
    quantities.push_back({"chi2_free_qop_self", "log_{10} #chi^{2} with q/p free (4 dof)", false,
                          {"self join", "same particle"}});
    quantities.push_back({"momentum_ratio", "|p_{b}| / |p_{a}| on the plane", false});
  }
  std::vector<double> chi2_free_qop_direction[2];  // same particle: forward, backward
  std::vector<double> chi2_direction[2];  // same particle: forward, backward

  auto find = [&](const TString& name) -> Quantity& {
    for (auto& quantity : quantities) {
      if (quantity.name == name) return quantity;
    }
    throw std::runtime_error(("unknown quantity " + name).Data());
  };

  // [class][ok]
  long n_pairs[2][2] = {{0, 0}, {0, 0}};
  long n_backward[2] = {0, 0};
  long n_self[2] = {0, 0};  // [ok]
  for (Long64_t i_entry = 0; i_entry < tree->GetEntries(); i_entry++) {
    tree->GetEntry(i_entry);
    if (self) {
      n_self[ok ? 1 : 0]++;
      if (ok && chi2 > 0) find("chi2_self").values[0].push_back(std::log10(chi2));
      if (has_free_qop && ok && chi2_free_qop > 0) find("chi2_free_qop_self").values[0].push_back(std::log10(chi2_free_qop));
      continue;
    }
    const int truth_class = (hash_a != 0 && hash_a == hash_b) ? 0 : 1;
    n_pairs[truth_class][ok ? 1 : 0]++;
    if (!ok) continue;
    n_backward[truth_class] += backward ? 1 : 0;
    if (chi2 > 0) {
      find("chi2").values[truth_class].push_back(std::log10(chi2));
      find("prob").values[truth_class].push_back(TMath::Prob(chi2, 5));
      if (truth_class == 0) {
        find("chi2_self").values[1].push_back(std::log10(chi2));
        chi2_direction[backward ? 1 : 0].push_back(chi2);
      }
    }
    if (has_free_qop) {
      if (chi2_free_qop > 0) {
        find("chi2_free_qop").values[truth_class].push_back(std::log10(chi2_free_qop));
        if (truth_class == 0) {
          find("chi2_free_qop_self").values[1].push_back(std::log10(chi2_free_qop));
          chi2_free_qop_direction[backward ? 1 : 0].push_back(chi2_free_qop);
        }
      }
      if (momentum_ratio > 0) find("momentum_ratio").values[truth_class].push_back(momentum_ratio);
    }
    const bool both_tof = (n_itof_a + n_otof_a > 0) && (n_itof_b + n_otof_b > 0);
    if (chi2_time > 0) {
      find(both_tof ? "chi2_time_tof" : "chi2_time_notof").values[truth_class].push_back(std::log10(chi2_time));
    }
    for (int i_par = 0; i_par < 6; i_par++) {
      const double sigma = std::hypot(sigma_a[i_par], sigma_b[i_par]);
      if (sigma > 0 && std::isfinite(delta[i_par])) {
        find(Form("pull_%d", i_par)).values[truth_class].push_back(delta[i_par] / sigma);
      }
    }
  }
  std::cout << "Track pairs (ok / total): same particle " << n_pairs[0][1] << " / " << n_pairs[0][0] + n_pairs[0][1]
            << ", other " << n_pairs[1][1] << " / " << n_pairs[1][0] + n_pairs[1][1] << std::endl;

  TString output_dir = Utils::aliceProjectsBase(onSTBC) + "acts-studies/low_pt_studies/figures/track_joins/"
                     + TString(input_dir).ReplaceAll("/", "_");
  gSystem->mkdir(output_dir, true);
  TFile* output = TFile::Open(output_dir + "/track_join_chi2.root", "RECREATE");

  TCanvas* canvas = new TCanvas("c_track_join_chi2", "Track joins: chi2 on the plane at the later track", 2400, 1500);
  const int n_rows = (quantities.size() + 1 + 3) / 4;  // + summary
  canvas->SetCanvasSize(2400, 500 * n_rows);
  canvas->Divide(4, n_rows);
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
    if (quantity.name == "prob") {
      low = 0;
      high = 1;
    }
    if (!(high > low)) high = low + 1;

    canvas->cd(i_quantity + 1);
    if (quantity.name.BeginsWith("pull")) gPad->SetLogy();
    TLegend* legend = new TLegend(0.5, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    double maximum = 0;
    TH1D* histograms[2];
    for (int i_class = 0; i_class < 2; i_class++) {
      TH1D* histogram = new TH1D(Form("h_%s_%d", quantity.name.Data(), i_class),
                                 Form(";%s;fraction of pairs", quantity.title.Data()), n_bins, low, high);
      for (double value : quantity.values[i_class]) histogram->Fill(value);
      histogram->SetLineColor(class_colors[i_class]);
      histogram->SetLineWidth(2);
      output->cd();
      histogram->Write();
      if (histogram->Integral(0, n_bins + 1) > 0) histogram->Scale(1. / histogram->Integral(0, n_bins + 1));
      maximum = std::max(maximum, histogram->GetMaximum());
      legend->AddEntry(histogram, Form("%s (%zu)", quantity.labels[i_class].Data(), quantity.values[i_class].size()), "l");
      histograms[i_class] = histogram;
    }
    for (int i_class = 0; i_class < 2; i_class++) {
      histograms[i_class]->SetMaximum(quantity.name.BeginsWith("pull") ? maximum * 3 : maximum * 1.3);
      if (quantity.name.BeginsWith("pull")) histograms[i_class]->SetMinimum(maximum * 1e-3);
      histograms[i_class]->Draw(i_class == 0 ? "hist" : "hist same");
    }
    // expected density of log10 chi2 with n dof, per bin: f(x) x ln(10) * bin width
    if (quantity.name == "chi2" || quantity.name == "chi2_self" || quantity.name.BeginsWith("chi2_free_qop")) {
      const double width = (high - low) / n_bins;
      const int dof = quantity.name.BeginsWith("chi2_free_qop") ? 4 : 5;
      TF1* expected = new TF1(Form("f_%s", quantity.name.Data()),
                              [width, dof](double* x, double*) {
                                const double chi2_value = std::pow(10., x[0]);
                                return ROOT::Math::chisquared_pdf(chi2_value, dof) * chi2_value * std::log(10.) * width;
                              },
                              low, high, 0);
      expected->SetLineColor(kBlack);
      expected->SetLineStyle(2);
      expected->SetNpx(500);
      expected->Draw("same");
      legend->AddEntry(expected, Form("#chi^{2}(%d)", dof), "l");
    }
    legend->Draw();
  }

  canvas->cd(4 * n_rows);
  TPaveText* summary = new TPaveText(0.05, 0.1, 0.95, 0.9, "NDC");
  summary->SetBorderSize(0);
  summary->SetFillStyle(0);
  summary->SetTextAlign(12);
  summary->AddText(Form("pairs: %ld", n_pairs[0][0] + n_pairs[0][1] + n_pairs[1][0] + n_pairs[1][1]));
  const char* class_names[2] = {"same particle", "other"};
  for (int truth_class = 0; truth_class < 2; truth_class++) {
    const long total = n_pairs[truth_class][0] + n_pairs[truth_class][1];
    summary->AddText(Form("%s: %ld, reached the plane %ld (%.0f%%), backward %ld", class_names[truth_class], total,
                          n_pairs[truth_class][1], total > 0 ? 100. * n_pairs[truth_class][1] / total : 0.,
                          n_backward[truth_class]));
  }
  for (int direction = 0; direction < 2; direction++) {
    summary->AddText(Form("same particle %s: %zu, median #chi^{2} %.3g", direction == 0 ? "forward" : "backward",
                          chi2_direction[direction].size(), quantile(chi2_direction[direction], 0.5)));
  }
  if (has_free_qop) {
    for (int direction = 0; direction < 2; direction++) {
      summary->AddText(Form("same particle %s: median #chi^{2} with q/p free %.3g",
                            direction == 0 ? "forward" : "backward",
                            quantile(chi2_free_qop_direction[direction], 0.5)));
    }
  }
  if (has_self) {
    summary->AddText(Form("self joins: %ld, reached the plane %ld, median #chi^{2} %.3g", n_self[0] + n_self[1],
                          n_self[1], std::pow(10., quantile(find("chi2_self").values[0], 0.5))));
  }
  summary->Draw();

  canvas->SaveAs(output_dir + "/track_join_chi2.pdf");
  output->cd();
  canvas->Write();
  output->Close();
  std::cout << "Wrote " << output_dir << "/track_join_chi2.{pdf,root}" << std::endl;
}
