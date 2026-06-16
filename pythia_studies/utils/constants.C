#pragma once

#include <TString.h>

namespace Constants {

constexpr double pi = 3.14159265358979323846;
constexpr double two_pi = 2.0 * pi;
constexpr double half_pi = 0.5 * pi;

const TString alice_base_local = "/home/justus/projects/alice";
const TString alice_base_stbc = "/data/alice/jrudolph/alice";

const TString getAliceBaseDir(bool onSTBC) {
  return onSTBC ? alice_base_stbc : alice_base_local;
}

}  // namespace Constants