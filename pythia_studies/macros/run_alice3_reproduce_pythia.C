#include <vector>
#include "alice3_reproduce.C"

void run_alice3_reproduce_pythia(
  const TString pathToFilesFromOutput,
  bool onSTBC=false, unsigned nFiles=1, bool include_muons=false) {
  std::vector<TString> input_dirs;
  for (unsigned i = 0; i < nFiles; i++) {
    input_dirs.push_back(Form("%s/seed_%u", pathToFilesFromOutput.Data(), i));
  }
  alice3_reproduce(onSTBC, input_dirs, input_dirs, 1.0, 0.5, include_muons);
}