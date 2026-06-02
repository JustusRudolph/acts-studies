#pragma once

#include <TString.h>

namespace Constants {

constexpr double pi = 3.14159265358979323846;
constexpr double two_pi = 2.0 * pi;
constexpr double half_pi = 0.5 * pi;

const TString actso2_base_local = "/home/justus/projects/alice/ACTSO2/";
const TString actso2_base_stbc = "/data/alice/jrudolph/alice/ACTSO2/";

const TString getACTSO2Base(bool onSTBC) {
  return onSTBC ? actso2_base_stbc : actso2_base_local;
}

}  // namespace Constants