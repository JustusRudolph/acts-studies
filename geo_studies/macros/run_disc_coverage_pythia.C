#include <vector>
#include "disc_coverage.C"

void run_disc_coverage_pythia(bool onSTBC=false,
                              float collision_rate=2400.,
                              unsigned nEvents=1000000,
                              bool runWithMeasurements=false,
                              unsigned nFiles=1,
                              float tolML=3.4,
                              float tolOT=3.4,
                              unsigned debugLevel=0) {
  std::vector<TString> input_dirs;
  for (unsigned i = 0; i < nFiles; i++) {
    input_dirs.push_back(Form("pythia_2500ev_pp_pileup0/seed_%u", i));
  }
  disc_coverage(input_dirs, onSTBC, collision_rate, nEvents, runWithMeasurements, tolML, tolOT, debugLevel);
}