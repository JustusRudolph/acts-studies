#pragma once

#include <cmath>
#include "constants.C"

namespace Cuts {

bool petalcut_phi(float phi, float epsilon=0.1) {
  // Cut out the petal regions in phi
  // Petal regions are at around -pi, -pi/2, 0, pi/2, pi
  float phi_abs = std::abs(phi);
  return (
          (   phi_abs < epsilon ) ||
          ( ( phi_abs > (Constants::half_pi - epsilon) ) && ( phi_abs < (Constants::half_pi + epsilon) ) ) ||
          ( ( phi_abs > (Constants::pi - epsilon) ) && ( phi_abs < (Constants::pi + epsilon) ) )
         );
}

bool eta_cut(float eta, float abs_eta_max) {
  return std::abs(eta) < abs_eta_max;
}

bool minHits_cut(unsigned nHits, unsigned min_nHits) {
  return nHits >= min_nHits;
}

bool primary_cut(unsigned generation) {
  return generation == 0;
}

bool alice3_default_cut(float eta, unsigned nHits, unsigned generation) {
  return eta_cut(eta, 1) &&
         minHits_cut(nHits, 7) &&
         primary_cut(generation);
}

}  // namespace Cuts