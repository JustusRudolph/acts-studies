#include <vector>
#include "disc_coverage_v2.C"

void run_disc_coverage_v2_pythia(TString pathToFilesFromOutput = "pythia_0pu_2500ev",
                                 unsigned lo_idx=0,
                                 unsigned hi_idx=0,
                                 bool onSTBC=false,
                                 float collision_rate=24000.,
                                 unsigned nEvents=1000000,
                                 bool runWithMeasurements=false,
                                 std::tuple<float, float> tolML={0, 2.5},
                                 std::tuple<float, float> tolOT={5, 3},
                                 unsigned debugLevel=0,
                                 bool gen3Geometry=true,
                                 TString tag="",  // names the histogram file
                                 Utils::AnalysisBase::Mode mode=Utils::AnalysisBase::kAuto) {
  std::vector<TString> input_dirs;
  for (unsigned i = lo_idx; i <= hi_idx; i++) {
    input_dirs.push_back(Form("%s/seed_%u", pathToFilesFromOutput.Data(), i));
  }
  disc_coverage_v2(input_dirs, onSTBC, collision_rate, nEvents, runWithMeasurements,
                   tolML, tolOT, debugLevel, gen3Geometry, tag, mode);
}
