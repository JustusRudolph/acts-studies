// Map of which (x, y) bins of the FT3 discs are covered by active silicon, built
// from the simulated geometry itself (o2sim_geometry.root). The O2 stave tables
// (FT3ModuleConstants.h) are not used on purpose: which module fill ends up in a
// geometry depends on the configKeyValues it was made with (e.g. whether
// FT3Base.useExactStavePlacement was set), the geometry file is what was simulated.
//
// There is one TH2C per disc, 1 where the bin centre lies on an FT3Sensor_Active
// volume and 0 elsewhere. Alongside it the staves (FT3_Stave volumes) of every disc
// are stored: x midpoint, x range of their active sensors and the z of those
// sensors. They give the stave IDs, so neither the number of staves nor their pitch
// has to be hard coded per geometry, and the z resolves which stave a hit is on
// where neighbouring staves overlap in x. Both are written to a ROOT file and read
// back from it on the next run, so the geometry only has to be walked once:
//
//   GeometryMapping::ActiveAreaMap map("histos/geo/active_area_map_tag.root",
//                                      ".../geom/o2sim_geometry.root");
//   map.setBinning(nXYBins, xyMax);     // must match the hit maps it is used with
//   if (!map.loadOrBuild(rebuild)) ...  // rebuild: ignore an existing map file
//   bool active = map.isActive(side, disc, x_mm, y_mm);
//   int stave = map.staveID(side, disc, x_mm, z_mm);
//
// Discs are indexed like in disc_coverage_v2.C: side 0 is backward (z < 0), side 1
// forward, and within a side the six discs go by increasing |z| (3 ML, then 3 OT).
#pragma once

#include <TFile.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TH2C.h>
#include <TString.h>
#include <TSystem.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace GeometryMapping {

class ActiveAreaMap {
 public:
  static constexpr int nSides = 2;  // 0 backward (z < 0), 1 forward (z > 0)
  static constexpr int nDiscs = 6;  // per side, by increasing |z|

  ActiveAreaMap(TString mapFilePath, TString geometryFilePath)
      : fMapFilePath(mapFilePath), fGeometryFilePath(geometryFilePath) {}
  ~ActiveAreaMap() { clear(); }
  // owns its histograms
  ActiveAreaMap(const ActiveAreaMap&) = delete;
  ActiveAreaMap& operator=(const ActiveAreaMap&) = delete;

  // square binning per disc, nBins in both x and y over [-xyMax, xyMax] (mm)
  void setBinning(const std::array<unsigned, nDiscs>& nBins,
                  const std::array<float, nDiscs>& xyMax_mm) {
    fNBins = nBins;
    fXYMax = xyMax_mm;
  }

  /*
   * Read the map file if it is there (and binned like setBinning() asks for),
   * otherwise build the map from the geometry and write it. rebuild skips the
   * reading, i.e. the map file is always overwritten.
   */
  bool loadOrBuild(bool rebuild) {
    if (!rebuild && read()) {
      printf("Read active area map from %s.\n", fMapFilePath.Data());
      return true;
    }
    if (!build()) return false;
    write();
    return true;
  }

  bool isActive(int side, int disc, float x_mm, float y_mm) const {
    const TH2C* hitmap = fMaps[side][disc];
    return hitmap && hitmap->GetBinContent(hitmap->FindFixBin(x_mm, y_mm)) > 0;
  }

  TH2C* map(int side, int disc) const { return fMaps[side][disc]; }

  // One stave of a disc (mm). The two pieces of a stave split by the beam pipe share
  // their x, so they are one stave here.
  struct Stave {
    double x;           // midpoint
    double xMin, xMax;  // x range covered by its active sensors
    double z;           // mean z of its active sensors, i.e. where its hits are
  };
  // the staves of a disc, sorted by x
  const std::vector<Stave>& staves(int side, int disc) const { return fStaves[side][disc]; }

  /*
   * Stave ID of a hit, counted like staveIdxToID() in FT3ModuleConstants.h: 1-indexed
   * from the middle outwards, negative for x < 0, no 0. The hit goes to the stave whose
   * active x range contains it; where neighbouring staves overlap in x (they sit at
   * different z, alternating front and back), to the one closest in z. A hit outside
   * every stave's x range (e.g. between sensors) goes to the closest midpoint.
   * inOverlap, if given, is set to whether the hit was in an overlap.
   * Returns 0 if the disc has no staves.
   */
  int staveID(int side, int disc, double x_mm, double z_mm, bool* inOverlap = nullptr) const {
    if (inOverlap) *inOverlap = false;
    const std::vector<Stave>& staves = fStaves[side][disc];
    if (staves.empty()) return 0;
    const int nStaves = staves.size();
    // closest midpoint
    auto it = std::lower_bound(staves.begin(), staves.end(), x_mm,
                               [](const Stave& s, double x) { return s.x < x; });
    if (it == staves.end() || (it != staves.begin() &&
                               ( x_mm - (it - 1)->x ) < ( it->x - x_mm ) ))
      --it;  // not at the first stave and closer to the previous x_mid
    int idx = it - staves.begin();
    // staves only overlap with their neighbours
    int best = idx, nContaining = 0;
    double bestDz = std::numeric_limits<double>::max();
    // check for overlap in neighbouring stave(s) - edge staves only have one neighbour
    for (int i = std::max(0, idx - 1); i <= std::min(nStaves - 1, idx + 1); i++) {
      if (x_mm < staves[i].xMin || x_mm > staves[i].xMax) continue;
      const double dz = std::abs(z_mm - staves[i].z);
      if (dz < bestDz) {
        best = i;
        bestDz = dz;
      }
      nContaining++;
    }
    if (inOverlap) *inOverlap = nContaining > 1;
    const int nNegative = nStavesNegativeX(side, disc);
    return best < nNegative ? best - nNegative : best - nNegative + 1;
  }

  // largest |stave ID| of a disc, e.g. to book a histogram from -N to +N
  // useful when the staves are not mirrored around x = 0
  int maxStaveID(int side, int disc) const {
    const int nNegative = nStavesNegativeX(side, disc);
    return std::max(nNegative, (int)fStaves[side][disc].size() - nNegative);
  }

  /*
   * Active area (mm^2) in each of nBins uniform r bins over [0, rMax_mm]. Map bins are
   * 0/1 flags and the hits inside a bin are not located any finer, so every active bin
   * puts its full area into the r bin of its centre.
   */
  std::vector<double> activeAreaVsR_mm2(int side, int disc, int nBins, double rMax_mm) const {
    std::vector<double> area(nBins, 0.);
    const TH2C* hitmap = fMaps[side][disc];
    if (!hitmap || nBins <= 0 || rMax_mm <= 0) return area;
    const double rBinWidth = rMax_mm / nBins;
    const TAxis* x_axis = hitmap->GetXaxis();
    const TAxis* y_axis = hitmap->GetYaxis();
    for (int x_bin = 1; x_bin <= x_axis->GetNbins(); x_bin++) {
      const double x = x_axis->GetBinCenter(x_bin), wx = x_axis->GetBinWidth(x_bin);
      for (int y_bin = 1; y_bin <= y_axis->GetNbins(); y_bin++) {
        if (hitmap->GetBinContent(x_bin, y_bin) <= 0) continue;
        const double y = y_axis->GetBinCenter(y_bin), wy = y_axis->GetBinWidth(y_bin);
        const int rBin = static_cast<int>(std::sqrt(x * x + y * y) / rBinWidth);
        if (rBin < nBins) area[rBin] += wx * wy;
      }
    }
    return area;
  }

 private:
  static TString mapName(int side, int disc) {
    return Form("activeMap_%s_disc_%d", side == 0 ? "bwd" : "fwd", disc);
  }
  // stored flat as x, xMin, xMax, z per stave
  static TString stavesName(int side, int disc) {
    return Form("staves_%s_disc_%d", side == 0 ? "bwd" : "fwd", disc);
  }
  static constexpr int nStaveFields = 4;

  int nStavesNegativeX(int side, int disc) const {
    const std::vector<Stave>& st = fStaves[side][disc];
    return std::count_if(st.begin(), st.end(), [](const Stave& s) { return s.x < 0; });
  }

  bool binningMatches(const TH2C* h, int disc) const {
    const double tol = 1e-3;  // mm
    return h->GetNbinsX() == (int)fNBins[disc] && h->GetNbinsY() == (int)fNBins[disc] &&
           std::abs(h->GetXaxis()->GetXmin() + fXYMax[disc]) < tol &&
           std::abs(h->GetXaxis()->GetXmax() - fXYMax[disc]) < tol &&
           std::abs(h->GetYaxis()->GetXmin() + fXYMax[disc]) < tol &&
           std::abs(h->GetYaxis()->GetXmax() - fXYMax[disc]) < tol;
  }

  void clear() {
    for (auto& side : fMaps)
      for (TH2C*& hitmap : side) {
        delete hitmap;
        hitmap = nullptr;
      }
    for (auto& side : fStaves)
      for (auto& st : side) st.clear();
  }

  bool read() {
    if (gSystem->AccessPathName(fMapFilePath)) return false;  // doesn't exist
    TFile* file = TFile::Open(fMapFilePath);
    if (!file || file->IsZombie()) return false;
    clear();
    bool ok = true;
    for (int side = 0; side < nSides && ok; side++) {
      for (int disc = 0; disc < nDiscs && ok; disc++) {
        TH2C* hitmap = dynamic_cast<TH2C*>(file->Get(mapName(side, disc)));
        TVectorD* st = dynamic_cast<TVectorD*>(file->Get(stavesName(side, disc)));
        if (!hitmap || !binningMatches(hitmap, disc) || !st || st->GetNrows() % nStaveFields != 0) {
          std::cout << "Active area map " << fMapFilePath << " is missing "
                    << mapName(side, disc) << "/" << stavesName(side, disc)
                    << " or binned differently - rebuilding it." << std::endl;
          delete hitmap;
          delete st;
          ok = false;
          break;
        }
        hitmap->SetDirectory(nullptr);
        fMaps[side][disc] = hitmap;
        for (int i = 0; i < st->GetNrows(); i += nStaveFields)
          fStaves[side][disc].push_back({(*st)[i], (*st)[i + 1], (*st)[i + 2], (*st)[i + 3]});
        delete st;
      }
    }
    file->Close();
    delete file;
    if (!ok) clear();
    return ok;
  }

  void write() const {
    gSystem->mkdir(gSystem->DirName(fMapFilePath), true);
    TFile* file = new TFile(fMapFilePath, "RECREATE");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not write the active area map to " << fMapFilePath
                << "." << std::endl;
      return;
    }
    for (int side = 0; side < nSides; side++) {
      for (int disc = 0; disc < nDiscs; disc++) {
        fMaps[side][disc]->Write();
        std::vector<double> flat;
        for (const Stave& st : fStaves[side][disc])
          flat.insert(flat.end(), {st.x, st.xMin, st.xMax, st.z});
        TVectorD(flat.size(), flat.data()).Write(stavesName(side, disc));
      }
    }
    file->Close();
    delete file;
    printf("Wrote active area map to %s.\n", fMapFilePath.Data());
  }

  /*
   * Identify the disc a sensor is on: the path down to the deepest FT3Layer volume
   * above it. NOTE: this has to be the placement (node path), not the volume name,
   * since the geometry reuses one layer volume for several discs (e.g. FT3Layer1_0
   * for all three forward ML discs), and node names don't contain the volume name.
   */
  static TString layerPlacement(const TGeoIterator& it) {
    int layerLevel = -1;
    for (int level = 1; level < it.GetLevel(); level++)
      if (TString(it.GetNode(level)->GetVolume()->GetName()).BeginsWith("FT3Layer"))
        layerLevel = level;
    if (layerLevel < 0) return "";
    TString placement;
    for (int level = 1; level <= layerLevel; level++)
      placement += TString("/") + it.GetNode(level)->GetName();
    // the volume name only to make the printout readable
    return placement + " (" + it.GetNode(layerLevel)->GetVolume()->GetName() + ")";
  }

  bool build() {
    if (gSystem->AccessPathName(fGeometryFilePath)) {
      std::cerr << "Error: Geometry file " << fGeometryFilePath << " not found." << std::endl;
      return false;
    }
    TGeoManager* geo = TGeoManager::Import(fGeometryFilePath);
    if (!geo) {
      std::cerr << "Error: Could not read a geometry from " << fGeometryFilePath << "."
                << std::endl;
      return false;
    }

    // xy footprint of every active sensor (mm), grouped by the disc layer it sits in
    struct Footprint { double xMin, xMax, yMin, yMax, z; };
    std::map<TString, std::vector<Footprint>> sensorsPerLayer;
    std::map<TString, double> zSumPerLayer;
    std::map<TString, std::vector<double>> staveXPerLayer;  // one entry per stave piece
    unsigned nNotAxisAligned = 0;
    TGeoIterator next(geo->GetTopVolume());
    while (TGeoNode* node = next()) {
      // match the volume, the node names of an imported GDML are just "<mother>_pv<n>"
      const TString volumeName = node->GetVolume()->GetName();
      const bool isSensor = volumeName.BeginsWith("FT3Sensor_Active");
      const bool isStave = volumeName.BeginsWith("FT3_Stave");
      if (!isSensor && !isStave) continue;
      const TString layer = layerPlacement(next);
      if (layer.IsNull()) continue;

      if (isStave) {  // only its position matters: where its local origin ends up
        const double local[3] = {0, 0, 0};
        double global[3];
        next.GetCurrentMatrix()->LocalToMaster(local, global);
        staveXPerLayer[layer].push_back(global[0] * 10);  // cm -> mm
        continue;
      }
      const TGeoBBox* box = dynamic_cast<const TGeoBBox*>(node->GetVolume()->GetShape());
      if (!box) continue;

      // the global bounding box of the 8 corners
      const TGeoMatrix* matrix = next.GetCurrentMatrix();
      const double half[3] = {box->GetDX(), box->GetDY(), box->GetDZ()};
      const double* origin = box->GetOrigin();
      double lo[3], hi[3];
      std::fill(lo, lo + 3, std::numeric_limits<double>::max());
      std::fill(hi, hi + 3, std::numeric_limits<double>::lowest());
      for (int corner = 0; corner < 8; corner++) {
        double local[3], global[3];
        for (int k = 0; k < 3; k++)
          local[k] = origin[k] + ((corner >> k) & 1 ? half[k] : -half[k]);
        matrix->LocalToMaster(local, global);
        for (int k = 0; k < 3; k++) {
          lo[k] = std::min(lo[k], global[k]);
          hi[k] = std::max(hi[k], global[k]);
        }
      }
      // the bounding box is only the footprint if the sensor is not turned in (x, y):
      // then its area equals that of the two large faces of the box
      std::array<double, 3> sorted = {half[0], half[1], half[2]};
      std::sort(sorted.begin(), sorted.end());
      const double faceArea = 4 * sorted[1] * sorted[2];
      if (std::abs((hi[0] - lo[0]) * (hi[1] - lo[1]) - faceArea) > 0.01 * faceArea)
        nNotAxisAligned++;

      // TGeo is in cm, the hits in mm
      sensorsPerLayer[layer].push_back(
        {lo[0] * 10, hi[0] * 10, lo[1] * 10, hi[1] * 10, 0.5 * (lo[2] + hi[2]) * 10});
      zSumPerLayer[layer] += 0.5 * (lo[2] + hi[2]) * 10;
    }
    delete geo;  // also resets gGeoManager

    if (nNotAxisAligned > 0)
      std::cerr << "Warning: " << nNotAxisAligned << " active sensors are turned in (x, y), "
                << "their footprint is overestimated by the bounding box." << std::endl;

    // order the layers per side by |z|
    std::array<std::vector<std::pair<double, TString>>, nSides> layersBySide;
    for (const auto& [layer, sensors] : sensorsPerLayer) {
      const double meanZ = zSumPerLayer[layer] / sensors.size();
      layersBySide[meanZ < 0 ? 0 : 1].push_back({std::abs(meanZ), layer});
    }
    for (int side = 0; side < nSides; side++) {
      if (layersBySide[side].size() != nDiscs) {
        std::cerr << "Error: Expected " << nDiscs << " FT3 layers with active sensors on the "
                  << (side == 0 ? "backward" : "forward") << " side of " << fGeometryFilePath
                  << ", found " << layersBySide[side].size() << "." << std::endl;
        return false;
      }
      std::sort(layersBySide[side].begin(), layersBySide[side].end());
    }

    // rasterise: a bin is active if its centre is on a sensor
    clear();
    for (int side = 0; side < nSides; side++) {
      for (int disc = 0; disc < nDiscs; disc++) {
        const auto& [absZ, layer] = layersBySide[side][disc];
        TH2C* h = new TH2C(mapName(side, disc),
                           Form("Active area %s (%s, |z| = %.1f mm);x (mm);y (mm)",
                                mapName(side, disc).Data(), layer.Data(), absZ),
                           fNBins[disc], -fXYMax[disc], fXYMax[disc],
                           fNBins[disc], -fXYMax[disc], fXYMax[disc]);
        h->SetDirectory(nullptr);
        h->SetStats(false);
        const TAxis* ax = h->GetXaxis();
        const TAxis* ay = h->GetYaxis();
        const int n = fNBins[disc];
        for (const Footprint& f : sensorsPerLayer[layer]) {
          const int bx0 = std::max(1, ax->FindFixBin(f.xMin));
          const int bx1 = std::min(n, ax->FindFixBin(f.xMax));
          const int by0 = std::max(1, ay->FindFixBin(f.yMin));
          const int by1 = std::min(n, ay->FindFixBin(f.yMax));
          for (int bx = bx0; bx <= bx1; bx++) {
            const double x = ax->GetBinCenter(bx);
            if (x < f.xMin || x > f.xMax) continue;
            for (int by = by0; by <= by1; by++) {
              const double y = ay->GetBinCenter(by);
              if (y < f.yMin || y > f.yMax) continue;
              h->SetBinContent(bx, by, 1);
            }
          }
        }
        fMaps[side][disc] = h;

        // the pieces of a split stave sit at the same x: keep one midpoint per stave
        std::vector<double> xs = staveXPerLayer[layer];
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end(),
                             [](double a, double b) { return std::abs(a - b) < 0.1; }),  // mm
                 xs.end());
        // every sensor belongs to the stave with the closest midpoint (its centre is well
        // inside the stave, so unlike a hit near the edge this is unambiguous), which
        // gives the active x range and z of each stave
        std::vector<Stave> st;
        for (double x : xs)
          st.push_back({x, std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::lowest(), 0.});
        std::vector<unsigned> nSensors(st.size(), 0);
        for (const Footprint& f : sensorsPerLayer[layer]) {
          if (st.empty()) break;
          const double xc = 0.5 * (f.xMin + f.xMax);
          auto closest = std::min_element(st.begin(), st.end(), [xc](const Stave& a, const Stave& b) {
            return std::abs(a.x - xc) < std::abs(b.x - xc);
          });
          closest->xMin = std::min(closest->xMin, f.xMin);
          closest->xMax = std::max(closest->xMax, f.xMax);
          closest->z += f.z;
          nSensors[closest - st.begin()]++;
        }
        unsigned nOverlaps = 0;
        double maxOverlap = 0;
        for (size_t i = 0; i < st.size(); i++) {
          if (nSensors[i] == 0) {  // no silicon: an empty x range, it can never be hit
            std::cerr << "Warning: stave at x = " << st[i].x << " mm of " << layer
                      << " has no active sensors." << std::endl;
            st[i].xMin = st[i].xMax = st[i].x;
            continue;
          }
          st[i].z /= nSensors[i];
          if (i > 0 && nSensors[i - 1] > 0 && st[i - 1].xMax > st[i].xMin) {
            nOverlaps++;
            maxOverlap = std::max(maxOverlap, st[i - 1].xMax - st[i].xMin);
          }
        }
        fStaves[side][disc] = st;

        printf("Active area map %s: %s at |z| = %.1f mm, %zu sensors, %.1f cm^2 active, "
               "%zu staves (%zu pieces), %u overlaps in x up to %.2f mm\n",
               mapName(side, disc).Data(), layer.Data(), absZ, sensorsPerLayer[layer].size(),
               h->Integral("width") / 100., st.size(), staveXPerLayer[layer].size(),
               nOverlaps, maxOverlap);
      }
    }
    return true;
  }

  TString fMapFilePath;
  TString fGeometryFilePath;
  std::array<unsigned, nDiscs> fNBins{};
  std::array<float, nDiscs> fXYMax{};
  std::array<std::array<TH2C*, nDiscs>, nSides> fMaps{};
  std::array<std::array<std::vector<Stave>, nDiscs>, nSides> fStaves;
};

}  // namespace GeometryMapping
