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

}  // namespace Cuts