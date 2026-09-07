// Base class for analysis macros. A macro derives from it, books its
// histograms in book(), fills them in fill() and draws them in plot(). The
// base takes care of writing every booked histogram out and reading them back
// in again, so that the (slow) loop over the data only has to be run once and
// plots can be iterated on quickly.
//
// It knows nothing about where the data comes from, opening it is part of the
// fill() implementation (see MCDataAnalysis.h for the LHCb MC checking files).
// A macro looks like this:
//
//   class MyAnalysis : public Utils::AnalysisBase {
//    public:
//     using Utils::AnalysisBase::AnalysisBase;
//     void book() override { bindHistogram<TH1D>("eta", "#eta;#eta;Counts", 50, -5., 6.); }
//     void fill() override { /* open the data and loop over it */ }
//     void plot() override { getHistFromFile<TH1D>("eta")->Draw(); }
//   };
//
//   void my_analysis(Utils::AnalysisBase::Mode mode=Utils::AnalysisBase::kAuto) {
//     MyAnalysis analysis("hists/my_hists.root");
//     analysis.run(mode);
//   }
//
// The first call fills the histogram file, every call after that reads it back
// and only redoes the plots. Pass kRefill to run over the data again:
//
//   root -l -q 'my_analysis.C(Utils::AnalysisBase::kRefill)'
#pragma once

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace Utils {

class AnalysisBase {
 public:
  // kAuto:     read the histogram file if it exists, otherwise fill it, then plot
  // kRefill:   always run over the data again, then plot
  // kFillOnly: run over the data and write the histograms, don't plot
  // kPlotOnly: only read the histogram file and plot (fails if it isn't there)
  enum Mode { kAuto, kRefill, kFillOnly, kPlotOnly };

  // histFilePath is where the booked histograms are written to and read back from.
  AnalysisBase(TString histFilePath) : fHistFilePath(histFilePath) {}

  virtual ~AnalysisBase() { clearObjects(); }

  // ---------------------------------------------------------------- macro API
  // Book every histogram here with bindHistogram<TH1D>("name", "title;x;y", ...).
  virtual void book() = 0;
  // Open whatever data is needed and fill the histograms booked above.
  virtual void fill() = 0;
  // Draw whatever is needed, the histograms are available through getHistFromFile<TH1D>().
  virtual void plot() = 0;

  /*
   * Drives the whole thing: get the histograms (from file or from the data)
   * and then plot them.
   */
  void run(Mode mode = kAuto) {
    applyStyle();
    bool needsFill = (mode == kRefill || mode == kFillOnly);
    if (!needsFill) {
      if (readHistsFromFile()) {
        printf("Read %zu histograms from %s.\n", fOrder.size(), fHistFilePath.Data());
        fBinding = true;  // let book() point the macro's members at what was read
        book();
        fBinding = false;
      } else if (mode == kPlotOnly) {
        std::cerr << "Error: No histogram file " << fHistFilePath
                  << ", run with kAuto or kRefill first." << std::endl;
        return;
      } else {
        needsFill = true;
      }
    }
    if (needsFill) {
      clearObjects();
      book();
      printf("Filling %zu histograms.\n", fOrder.size());
      fill();
      writeHistsToFile();
    }
    if (mode != kFillOnly) plot();
  }

  // ------------------------------------------------------- histogram registry
  /*
   * Bind a histogram (or any other TObject) to the name it is stored under:
   * creates it and keeps track of it, so that it is written out and read back
   * in without having to list it anywhere else.
   *
   * After the histograms have been read back from file, book() is called again
   * in binding mode, where this hands out the object that was read instead of
   * creating a new one. That way a macro that keeps its histograms in members
   * (rather than looking them up by name every time) can just write
   *   hEta = bindHistogram<TH1D>("eta", ...);
   * and have the member point at the right thing in both cases.
   */
  template <class T, class... Args>
  T* bindHistogram(const TString& name, Args&&... args) {
    if (fBinding) return getHistFromFile<T>(name);  // already read from the histogram file
    T* obj = new T(name.Data(), std::forward<Args>(args)...);
    adoptNewHistogram(obj);
    return obj;
  }

  // Same, for objects that can't be created with the templated version above.
  void adoptNewHistogram(TObject* obj) {
    if (!obj) return;
    detachHistogram(obj);
    // the stat box is off for everything unless setShowStats(true) was called
    // SetStats is only a function for TH1 derived objects, so check first if
    // the object is a TH1 before setting
    if (TH1* hist = dynamic_cast<TH1*>(obj)) hist->SetStats(fShowStats);
    TString name = obj->GetName();
    if (fObjects.count(name)) {
      std::cerr << "Warning: Overwriting already booked object " << name << "." << std::endl;
      delete fObjects[name];
    } else {
      fOrder.push_back(name);
    }
    fObjects[name] = obj;
  }

  // Get a booked/read histogram back, prints a warning if it isn't there.
  template <class T>
  T* getHistFromFile(const TString& name) const {
    auto it = fObjects.find(name);
    if (it == fObjects.end()) {
      std::cerr << "Warning: No object called " << name << " found." << std::endl;
      return nullptr;
    }
    T* obj = dynamic_cast<T*>(it->second);
    if (!obj)
      std::cerr << "Warning: " << name << " is a " << it->second->ClassName()
                << ", not the requested type." << std::endl;
    return obj;
  }

  // Write every booked object to the histogram file (defaults to the one the
  // analysis was constructed with).
  void writeHistsToFile(TString path = "") {
    if (path.IsNull()) path = fHistFilePath;
    gSystem->mkdir(gSystem->DirName(path), true);
    TFile* outFile = new TFile(path, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
      std::cerr << "Error: Could not write to " << path << "." << std::endl;
      return;
    }
    outFile->cd();
    for (const TString& name : fOrder) fObjects[name]->Write(name);
    outFile->Close();
    delete outFile;
    printf("Wrote %zu histograms to %s.\n", fOrder.size(), path.Data());
  }

  /*
   * Read everything back from the histogram file, whatever is in there. The
   * objects are detached from the file so they stay alive after closing it.
   */
  bool readHistsFromFile(TString path = "") {
    if (path.IsNull()) path = fHistFilePath;
    if (gSystem->AccessPathName(path)) return false;  // doesn't exist
    TFile* file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    clearObjects();
    TIter nextKey(file->GetListOfKeys());
    while (TKey* key = (TKey*) nextKey()) {
      TString name = key->GetName();
      if (fObjects.count(name)) continue;  // no duplicate names
      TObject* obj = file->Get(name);
      if (obj) adoptNewHistogram(obj);
    }
    file->Close();
    delete file;
    return !fOrder.empty();
  }

  // ------------------------------------------------------------------ helpers
  // Uniform bin edges, as they are needed for pretty much every histogram.
  static std::vector<float> binEdges(unsigned nBins, float low, float high) {
    std::vector<float> edges(nBins + 1);
    float step = (high - low) / nBins;
    for (unsigned i = 0; i <= nBins; i++) edges[i] = low + step * i;
    return edges;
  }

  // Save a canvas to path, creating the directory if it isn't there yet.
  void saveCanvas(TCanvas* canvas, TString path) {
    gSystem->mkdir(gSystem->DirName(path), true);
    canvas->SaveAs(path);
  }

  // Print progress every kStep entries of a loop over nTotal.
  void printProgress(Long64_t entry, Long64_t nTotal, const char* what,
                     Long64_t kStep = 1000000) const {
    if (entry && (entry % kStep == 0))
      printf("Finished with %d%% of %s.\n", (int) (100. * (double) entry / nTotal), what);
  }

  TString histFilePath() const { return fHistFilePath; }
  void setHistFilePath(TString path) { fHistFilePath = path; }

  /*
   * Whether histograms carry a stat box, off by default. Unlike the global
   * gStyle->SetOptStat() in applyStyle() this is a property of the histogram
   * itself, so it is written to the histogram file and survives being read
   * back. Call before book() (i.e. before run()) to affect everything.
   */
  void setShowStats(bool show) { fShowStats = show; }

 protected:
  // Batch mode and no stat box, override for a different style.
  virtual void applyStyle() {
    gROOT->SetBatch();  // so stuff isn't autoplotted
    gStyle->SetOptStat(0);  // remove the info box
    gStyle->SetAxisMaxDigits(3);
  }

 private:
  // Objects owned by us, not by whatever file happens to be open.
  static void detachHistogram(TObject* obj) {
    // TObject does not have a SetDirectory, so we must check the type of the object
    // before detaching it. NOTE: Add more here if we get more types with SetDirectory.
    if (TH1* hist = dynamic_cast<TH1*>(obj)) hist->SetDirectory(nullptr);
    else if (TEfficiency* eff = dynamic_cast<TEfficiency*>(obj)) eff->SetDirectory(nullptr);
    else {
      std::cerr << "Warning: Tried to detach an object of unsupported type: "
                << obj->ClassName() << std::endl;
    }
  }

  void clearObjects() {
    for (const TString& name : fOrder) delete fObjects[name];
    fObjects.clear();
    fOrder.clear();
  }

  TString fHistFilePath;
  bool fBinding{false};                    // book() looks up instead of creating
  bool fShowStats{false};                  // stat box on booked/read histograms
  std::vector<TString> fOrder;             // booked objects, in booking order
  std::map<TString, TObject*> fObjects;
};

}  // namespace Utils
