#include <vector>
#include "geo_efficiency_plots.C"

void run_geo_eff_pythia(bool onSTBC=false,
                              unsigned nFiles=1,
                              const unsigned nMinHits=7,
                              const float etaMin=-3,
                              const float etaMax=3) {
  std::vector<TString> input_dirs;
  for (unsigned i = 0; i < nFiles; i++) {
    input_dirs.push_back(Form("pythia_2500ev/seed_%u", i));
  }
  geo_efficiency_plots(input_dirs, onSTBC, nMinHits, etaMin, etaMax);
}