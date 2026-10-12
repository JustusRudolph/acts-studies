/*
 * BasicFunctions.h
 *
 * This file contains basic utility functions used throughout the analysis macros.
 */

#pragma once

#include <TString.h>
#include <TSystem.h>

#include <vector>

namespace Utils {
  // function for getting the path base of output material
  inline TString outputPathBase(bool onSTBC) {
    return onSTBC ? "/dcache/alice/jrudolph/alice/ACTSO2/output/"
                  : "/home/justus/projects/alice/ACTSO2/output/";
  }
  inline TString aliceProjectsBase(bool onSTBC) {
    return onSTBC ? "/data/alice/jrudolph/alice/"
                  : "/home/justus/projects/alice/";
  }
  // directory the figures of a run go into: the tag if given, else the sample the
  // input paths belong to (e.g. "pythia_0pu_2500ev" for "pythia_0pu_2500ev/seed_3")
  inline TString figureDirName(const TString& tag, const std::vector<TString>& pathBases) {
    if (!tag.IsNull()) return tag;
    TString sample = pathBases.empty() ? TString("unknown") : pathBases.front();
    if (sample.Contains("/")) sample = gSystem->DirName(sample);
    sample.ReplaceAll("/", "_");
    return sample;
  }
}  // namespace Utils
