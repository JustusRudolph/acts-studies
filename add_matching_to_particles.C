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
                         + "/particles_digitized_selected.root";
  const TString perfAmbiFile = base + pathToFilesFromOutput
                         + "/performance_finding_ambi.root";
  const TString perfSeedFile = base + pathToFilesFromOutput
                         + "/performance_seeding.root";
  const TString outFile  = base + pathToFilesFromOutput
                         + "/particles_matched.root";
  // --- Step 1: build lookup map from matchingdetails ---
  TFile* fPerf = TFile::Open(perfAmbiFile);
  if (!fPerf || fPerf->IsZombie()) {
    std::cerr << "Cannot open " << perfAmbiFile << std::endl;
    return;
  }
  TFile* fSeedPerf = TFile::Open(perfSeedFile);
  if (!fSeedPerf || fSeedPerf->IsZombie()) {
    std::cerr << "Cannot open " << perfSeedFile << std::endl;
    return;
  }

  std::cout << "Opened files: " << perfAmbiFile << " and " << perfSeedFile << std::endl;
  TTree* tPerf = (TTree*)fPerf->Get("matchingdetails");
  if (!tPerf) { std::cerr << "matchingdetails tree not found" << std::endl; return; }

  TTree* tSeedPerf = (TTree*)fSeedPerf->Get("matchingdetails");
  if (!tSeedPerf) {
    std::cerr << "matchingdetails tree not found in seeding file" << std::endl;
    return;
  }

  // particle is a scalar value in matchingdetails
  uint32_t perfParticle;
  uint32_t perfEventId;
  std::vector<uint32_t>* matchedIdxs     = nullptr;
  std::vector<double>*   matchedWeights  = nullptr;
  std::vector<uint32_t>* fakeIdxs        = nullptr;
  std::vector<double>*   fakeWeights     = nullptr;

  uint32_t seedPerfParticle;
  uint32_t seedPerfEventId;
  std::vector<uint32_t>* seedMatchedIdxs     = nullptr;
  std::vector<double>*   seedMatchedWeights  = nullptr;
  std::vector<uint32_t>* seedFakeIdxs        = nullptr;
  std::vector<double>*   seedFakeWeights     = nullptr;

  // set other branch addresses later to avoid loading vectors for this step
  tPerf->SetBranchAddress("event_nr", &perfEventId);
  tPerf->SetBranchAddress("particle_id_particle", &perfParticle);
  tSeedPerf->SetBranchAddress("event_nr", &seedPerfEventId);
  tSeedPerf->SetBranchAddress("particle_id_particle", &seedPerfParticle);

  // Nth entry in the vector corresponds to event N
  // map is (particle_id) ->  idx of matched & fake tracks in global matchingdetails tree
  using MatchMap = std::unordered_map<unsigned, std::map<unsigned, unsigned>>;
  MatchMap matchMap(tPerf->GetEntries());
  MatchMap seedMatchMap(tSeedPerf->GetEntries());

  std::cout << "Loading matchingdetails from " << perfAmbiFile <<
            " with " << tPerf->GetEntries() << " entries." << std::endl;

  std::cout << "matchMap size: " << matchMap.size() << std::endl;

  for (Long64_t i = 0; i < tPerf->GetEntries(); i++) {
    tPerf->GetEntry(i);
    matchMap[perfEventId] = {};
    matchMap[perfEventId][perfParticle] = i;
  }
  std::cout << "Loaded " << matchMap.size()
            << " entries from ambi matchingdetails" << std::endl;

  std::cout << "Loading matchingdetails from " << perfSeedFile <<
            " with " << tSeedPerf->GetEntries() << " entries." << std::endl;
  for (Long64_t i = 0; i < tSeedPerf->GetEntries(); i++) {
    tSeedPerf->GetEntry(i);
    seedMatchMap[seedPerfEventId] = {};
    seedMatchMap[seedPerfEventId][seedPerfParticle] = i;
  }
  std::cout << "Loaded " << seedMatchMap.size()
            << " entries from seeding matchingdetails" << std::endl;

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
  tSeedPerf->SetBranchAddress("matched_track_idxs",     &seedMatchedIdxs);
  tSeedPerf->SetBranchAddress("matched_track_weights",  &seedMatchedWeights);
  tSeedPerf->SetBranchAddress("fake_track_idxs",        &seedFakeIdxs);
  tSeedPerf->SetBranchAddress("fake_track_weights",     &seedFakeWeights);
  std::vector<std::vector<unsigned>> outMatchedIdxs;
  std::vector<std::vector<double>>   outMatchedWeights;
  std::vector<std::vector<unsigned>> outFakeIdxs;
  std::vector<std::vector<double>>   outFakeWeights;
  std::vector<std::vector<unsigned>> outSeedMatchedIdxs;
  std::vector<std::vector<double>>   outSeedMatchedWeights;
  std::vector<std::vector<unsigned>> outSeedFakeIdxs;
  std::vector<std::vector<double>>   outSeedFakeWeights;
  tOut->Branch("matched_track_idxs",    &outMatchedIdxs);
  tOut->Branch("matched_track_weights", &outMatchedWeights);
  tOut->Branch("fake_track_idxs",       &outFakeIdxs);
  tOut->Branch("fake_track_weights",    &outFakeWeights);
  tOut->Branch("seed_matched_track_idxs",    &outSeedMatchedIdxs);
  tOut->Branch("seed_matched_track_weights", &outSeedMatchedWeights);
  tOut->Branch("seed_fake_track_idxs",       &outSeedFakeIdxs);
  tOut->Branch("seed_fake_track_weights",    &outSeedFakeWeights);

  Long64_t nSim   = tSim->GetEntries();
  Long64_t nFound = 0;

  for (Long64_t i = 0; i < nSim; i++) {
    tSim->GetEntry(i);
    if (i % 1000 == 0 && i) {
      std::cout << "Processed " << i << " / " << nSim
                << " events. Found matches for " << nFound
                << " particles." << std::endl;
    }
    
    if (!simParticle->empty()) {  // There are particles in this event
      for (unsigned i_part = 0; i_part < simParticle->size(); i_part++) {
        unsigned simPID = (*simParticle)[i_part];
        if (checkPrimaryOnly && (*simGeneration)[i_part] != 0) {
          // If we're only checking primaries, skip particle if generation is not 0
          outMatchedIdxs.push_back({});
          outMatchedWeights.push_back({});
          outFakeIdxs.push_back({});
          outFakeWeights.push_back({});
          outSeedMatchedIdxs.push_back({});
          outSeedMatchedWeights.push_back({});
          outSeedFakeIdxs.push_back({});
          outSeedFakeWeights.push_back({});
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
        // now add seeding eff
        auto itSeedEvent = seedMatchMap.find(simEventId);
        if (itSeedEvent != seedMatchMap.end()) { // at least one particle has match
          auto itSeedPID = itSeedEvent->second.find(simPID);
          if (itSeedPID != itSeedEvent->second.end()) {
            unsigned seedPerfIdx = itSeedPID->second;
            tSeedPerf->GetEntry(seedPerfIdx);
            outSeedMatchedIdxs.push_back(*seedMatchedIdxs);
            outSeedMatchedWeights.push_back(*seedMatchedWeights);
            outSeedFakeIdxs.push_back(*seedFakeIdxs);
            outSeedFakeWeights.push_back(*seedFakeWeights);
          } else { // no match for this particle id
            outSeedMatchedIdxs.push_back({});
            outSeedMatchedWeights.push_back({});
            outSeedFakeIdxs.push_back({});
            outSeedFakeWeights.push_back({});
          }
        }  // else: no match for this simulated particle: fill empty vectors
        else {
          outSeedMatchedIdxs.push_back({});
          outSeedMatchedWeights.push_back({});
          outSeedFakeIdxs.push_back({});
          outSeedFakeWeights.push_back({});
        }
      }  // loop over simulated particles in event
    }
    tOut->Fill();

    outMatchedIdxs.clear();
    outMatchedWeights.clear();
    outFakeIdxs.clear();
    outFakeWeights.clear();
    outSeedMatchedIdxs.clear();
    outSeedMatchedWeights.clear();
    outSeedFakeIdxs.clear();
    outSeedFakeWeights.clear();
  }

  fOut->Write();

  fPerf->Close();
  fSeedPerf->Close();
  fOut->Close();
  fSim->Close();

  std::cout << "Matched " << nFound << " / " << nSim << " particles." << std::endl;
  std::cout << "Saved: " << outFile << std::endl;
}
