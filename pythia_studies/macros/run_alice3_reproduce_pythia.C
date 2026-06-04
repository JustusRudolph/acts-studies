#include <vector>
#include "alice3_reproduce.C"

void run_alice3_reproduce_pythia(bool onSTBC=false, unsigned nFiles=1) {
  std::vector<TString> input_dirs;
  for (unsigned i = 0; i < nFiles; i++) {
    input_dirs.push_back(Form("pythia_5k/seed_%u", i));
  }
  alice3_reproduce(onSTBC, input_dirs, input_dirs);
}