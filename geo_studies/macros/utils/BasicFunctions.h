/*
 * BasicFunctions.h
 *
 * This file contains basic utility functions used throughout the analysis macros.
 */

#pragma once

#include <TString.h>

namespace Utils {
  // function for getting the path base of output material
  inline TString outputPathBase(bool onSTBC) {
    return onSTBC ? "/dchache/alice/jrudolph/alice/ACTSO2/output/"
                  : "/home/justus/projects/alice/ACTSO2/output/";
  }
  inline TString aliceProjectsBase(bool onSTBC) {
    return onSTBC ? "/data/alice/jrudolph/alice/"
                  : "/home/justus/projects/alice/";
  }
}  // namespace Utils