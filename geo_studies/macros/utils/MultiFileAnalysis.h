// Base class for analyses that run over a directory of directories: one sub
// directory per job/seed, each holding several ROOT files (hits.root,
// measurements.root, ...). It adds the path handling and the "was this file
// written completely" gate to AnalysisBase, so that a job which was killed
// half way through never contributes an unnormalisable fraction of its events.
//
//   class MyAnalysis : public Utils::MultiFileAnalysis {
//     void fill() override {
//       for (const TString& subDir : subDirs()) {
//         TFile* hits = openFile(subDir, "hits.root");
//         if (!hits) continue;  // missing, unreadable or incomplete
//         TTree* tree = getTree(hits, "hits");
//         ...
//         hits->Close();
//       }
//     }
//   };
#pragma once

#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "AnalysisBase.h"

namespace Utils {

class MultiFileAnalysis : public AnalysisBase {
 public:
  /*
   * inputBase is the directory holding the sub directories, subDirs are the
   * ones to run over (leave empty and call discoverSubDirs() to take all of
   * them). Sub directories may themselves contain a "/", e.g. "pythia/seed_0".
   */
  MultiFileAnalysis(TString histFilePath, TString inputBase = "",
                    std::vector<TString> subDirs = {})
      : AnalysisBase(histFilePath), fInputBase(withSlash(inputBase)),
        fSubDirs(subDirs) {}

  void setInputBase(TString base) { fInputBase = withSlash(base); }
  void setSubDirs(std::vector<TString> subDirs) { fSubDirs = subDirs; }
  const std::vector<TString>& subDirs() const { return fSubDirs; }
  TString inputBase() const { return fInputBase; }

  /*
   * Take every sub directory of the input base, optionally only those whose
   * name starts with prefix, sorted by name (so seed_10 comes before seed_2).
   * Returns how many were found.
   */
  unsigned discoverSubDirs(TString prefix = "") {
    fSubDirs.clear();
    TSystemDirectory base(fInputBase, fInputBase);
    TList* entries = base.GetListOfFiles();
    if (!entries) {
      std::cerr << "Error: Could not list " << fInputBase << "." << std::endl;
      return 0;
    }
    TIter nextEntry(entries);
    while (TSystemFile* entry = (TSystemFile*) nextEntry()) {
      TString name = entry->GetName();
      if (!entry->IsDirectory() || name == "." || name == "..") continue;
      if (!prefix.IsNull() && !name.BeginsWith(prefix)) continue;
      fSubDirs.push_back(name);
    }
    delete entries;
    std::sort(fSubDirs.begin(), fSubDirs.end());
    return fSubDirs.size();
  }

  TString dirPath(const TString& subDir) const { return fInputBase + subDir + "/"; }
  TString filePath(const TString& subDir, const TString& fileName) const {
    return dirPath(subDir) + fileName;
  }

  /*
   * Open one of the files in a sub directory. Returns nullptr, with a message,
   * if it is missing, unreadable, or was not closed properly by the writer.
   * The caller owns the file and should Close() it when done.
   */
  TFile* openFile(const TString& subDir, const TString& fileName,
                  bool requireComplete = true) {
    TString path = filePath(subDir, fileName);
    // reject before opening: an incomplete file has an unknowable event count,
    // so it is dropped whole rather than contributing a partial sample
    if (requireComplete && !fileClosedProperly(path)) {
      std::cerr << "File not closed properly: " << path << " - skipping" << std::endl;
      return nullptr;
    }
    TFile* file = TFile::Open(path);
    if (!file || file->IsZombie()) {
      std::cerr << "Cannot open " << path << " - skipping" << std::endl;
      return nullptr;
    }
    // backstop: the header check is a heuristic, so drop anything ROOT still
    // had to recover
    if (requireComplete && file->TestBit(TFile::kRecovered)) {
      std::cerr << "File was recovered by ROOT: " << path << " - skipping" << std::endl;
      file->Close();
      delete file;
      return nullptr;
    }
    return file;
  }

  static TTree* getTree(TFile* file, const char* treeName) {
    if (!file) return nullptr;
    TTree* tree = (TTree*) file->Get(treeName);
    if (!tree)
      std::cerr << treeName << " tree not found in " << file->GetName() << std::endl;
    return tree;
  }

  /*
   * Read the ROOT file header directly to decide whether the writer closed the
   * file. A job that was killed mid-write never gets its free-segment record,
   * so fSeekFree stays 0. Checking this before TFile::Open means ROOT never
   * runs its (slow) key recovery, and we never read half-written events.
   */
  static bool fileClosedProperly(const TString& path) {
    std::ifstream in(path.Data(), std::ios::binary);
    if (!in) return false;
    char magic[4];
    in.read(magic, 4);
    if (in.gcount() != 4 || std::strncmp(magic, "root", 4) != 0) return false;

    auto readBE = [&in](std::streamoff off, int nbytes) -> uint64_t {
      in.seekg(off, std::ios::beg);
      unsigned char b[8];
      in.read(reinterpret_cast<char*>(b), nbytes);
      if (in.gcount() != nbytes) return 0;
      uint64_t v = 0;
      for (int i = 0; i < nbytes; i++) v = (v << 8) | b[i];  // ROOT headers are big-endian
      return v;
    };

    // header layout: fVersion at 4, fBEGIN at 8, fEND at 12, then fSeekFree.
    // fEND/fSeekFree take 8 bytes instead of 4 once the file is in large-file format.
    const bool largeFile = readBE(4, 4) > 1000000;  // fVersion
    const uint64_t seekFree   = largeFile ? readBE(20, 8) : readBE(16, 4);
    const uint64_t nbytesFree = largeFile ? readBE(28, 4) : readBE(20, 4);
    return seekFree > 0 && nbytesFree > 0;
  }

 private:
  static TString withSlash(TString dir) {
    if (!dir.IsNull() && !dir.EndsWith("/")) dir += "/";
    return dir;
  }

  TString fInputBase;
  std::vector<TString> fSubDirs;
};

}  // namespace Utils
