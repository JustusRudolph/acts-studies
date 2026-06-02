#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

void add_matching_to_particles(
  const TString pathToFilesFromOutput = "geo_staves_pi_1GeV_eta14-20",
  const bool onSTBC = false,
  const bool checkPrimaryOnly=true)
{

  TString base = "/home/justus/projects/alice/ACTSO2/output/";
  if (onSTBC) {
    base = "/data/alice/jrudolph/alice/ACTSO2/output/";
  }
  const TString simFile  = base + pathToFilesFromOutput
                         + "/particles_simulation.root";
  const TString perfFile = base + pathToFilesFromOutput
                         + "/performance_finding_ambi.root";
  const TString outFile  = base + pathToFilesFromOutput
                         + "/particles_simulation_matched.root";
  // --- Step 1: build lookup map from matchingdetails ---
  TFile* fPerf = TFile::Open(perfFile);
  if (!fPerf || fPerf->IsZombie()) { std::cerr << "Cannot open " << perfFile << std::endl; return; }

  TTree* tPerf = (TTree*)fPerf->Get("matchingdetails");
  if (!tPerf) { std::cerr << "matchingdetails tree not found" << std::endl; return; }

  // particle is a scalar value in matchingdetails
  uint32_t perfParticle;
  uint32_t perfEventId;
  std::vector<uint32_t>* matchedIdxs     = nullptr;
  std::vector<double>*   matchedWeights  = nullptr;
  std::vector<uint32_t>* fakeIdxs        = nullptr;
  std::vector<double>*   fakeWeights     = nullptr;

  // set other branch addresses later to avoid loading vectors for this step
  tPerf->SetBranchAddress("event_nr", &perfEventId);
  tPerf->SetBranchAddress("particle_id_particle", &perfParticle);

  // Nth entry in the vector corresponds to event N
  // map is (particle_id) ->  idx of matched & fake tracks in global matchingdetails tree
  std::unordered_map<unsigned, std::map<unsigned, unsigned>> matchMap(tPerf->GetEntries());

  std::cout << "Loading matchingdetails from " << perfFile <<
            " with " << tPerf->GetEntries() << " entries." << std::endl;

  std::cout << "matchMap size: " << matchMap.size() << std::endl;

  for (Long64_t i = 0; i < tPerf->GetEntries(); i++) {
    tPerf->GetEntry(i);
    // if (i > 9300) {
    //   std::cout << "Processing entry " << i << " with event_id " << perfEventId
    //             << " and particle_id " << perfParticle << std::endl;
    // }
    matchMap[perfEventId] = {};
    matchMap[perfEventId][perfParticle] = i;
  }

  std::cout << "Loaded " << matchMap.size() << " entries from matchingdetails" << std::endl;

  // --- Step 2: copy particles tree and add new branches ---
  TFile* fSim = TFile::Open(simFile);
  if (!fSim || fSim->IsZombie()) { std::cerr << "Cannot open " << simFile << std::endl; return; }

  TTree* tSim = (TTree*)fSim->Get("particles");
  if (!tSim) { std::cerr << "particles tree not found in " << simFile << std::endl; return; }

  TFile* fOut = TFile::Open(outFile, "RECREATE");
  TTree* tOut = tSim->CloneTree(0);

  // particle is a vector in particles_simulation, per event
  std::vector<unsigned>* simParticle = nullptr;
  std::vector<unsigned>* simGeneration = nullptr;
  unsigned simEventId;

  tSim->SetBranchAddress("event_id", &simEventId);
  tSim->SetBranchAddress("particle", &simParticle);
  tSim->SetBranchAddress("generation", &simGeneration);

  tPerf->SetBranchAddress("matched_track_idxs",     &matchedIdxs);
  tPerf->SetBranchAddress("matched_track_weights",  &matchedWeights);
  tPerf->SetBranchAddress("fake_track_idxs",        &fakeIdxs);
  tPerf->SetBranchAddress("fake_track_weights",     &fakeWeights);
  std::vector<std::vector<unsigned>> outMatchedIdxs;
  std::vector<std::vector<double>>   outMatchedWeights;
  std::vector<std::vector<unsigned>> outFakeIdxs;
  std::vector<std::vector<double>>   outFakeWeights;
  tOut->Branch("matched_track_idxs",    &outMatchedIdxs);
  tOut->Branch("matched_track_weights", &outMatchedWeights);
  tOut->Branch("fake_track_idxs",       &outFakeIdxs);
  tOut->Branch("fake_track_weights",    &outFakeWeights);

  Long64_t nSim   = tSim->GetEntries();
  Long64_t nFound = 0;

  for (Long64_t i = 0; i < nSim; i++) {
    tSim->GetEntry(i);
    
    if (!simParticle->empty()) {  // There are particles in this event
      for (unsigned i_part = 0; i_part < simParticle->size(); i_part++) {
        unsigned simPID = (*simParticle)[i_part];
        if (checkPrimaryOnly && (*simGeneration)[i_part] != 0) {
          // If we're only checking primaries, skip particle if generation is not 0
          outMatchedIdxs.push_back({});
          outMatchedWeights.push_back({});
          outFakeIdxs.push_back({});
          outFakeWeights.push_back({});
          continue;
        }
        auto itEvent = matchMap.find(simEventId);
        if (itEvent != matchMap.end()) { // at least one particle has match
          auto itPID = itEvent->second.find(simPID);
          if (itPID != itEvent->second.end()) {
            unsigned perfIdx = itPID->second;
            tPerf->GetEntry(perfIdx);
            outMatchedIdxs.push_back(*matchedIdxs);
            outMatchedWeights.push_back(*matchedWeights);
            outFakeIdxs.push_back(*fakeIdxs);
            outFakeWeights.push_back(*fakeWeights);
            // std::cout << "Primary matched index for ev " << simEventId << " particle "
            //           << simPID << ": " << (*matchedIdxs)[0]
            //           << " with weight " << (*matchedWeights)[0] << std::endl;
            nFound++;
          } else { // no match for this particle id
            outMatchedIdxs.push_back({});
            outMatchedWeights.push_back({});
            outFakeIdxs.push_back({});
            outFakeWeights.push_back({});
          }
        }  // else: no match for this simulated particle: fill empty vectors
        else {
          outMatchedIdxs.push_back({});
          outMatchedWeights.push_back({});
          outFakeIdxs.push_back({});
          outFakeWeights.push_back({});
        }
      }  // loop over simulated particles in event
    }
    tOut->Fill();

    outMatchedIdxs.clear();
    outMatchedWeights.clear();
    outFakeIdxs.clear();
    outFakeWeights.clear();
  }

  fOut->Write();

  fPerf->Close();
  fOut->Close();
  fSim->Close();

  std::cout << "Matched " << nFound << " / " << nSim << " particles." << std::endl;
  std::cout << "Saved: " << outFile << std::endl;
}
