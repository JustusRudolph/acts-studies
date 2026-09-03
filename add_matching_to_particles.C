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
  const bool useIterativeTracking=true,  // changes name of perf files
  const bool checkPrimaryOnly=true,
  const bool useDigitisedParticles=false)
{

  TString base = "/home/justus/projects/alice/ACTSO2/output/";
  if (onSTBC) {
    base = "/data/alice/jrudolph/alice/ACTSO2/output/";
  }
  const TString simFile  = base + pathToFilesFromOutput
                         + (useDigitisedParticles ?
                              "/particles_digitized_selected.root" :
                              "/particles_simulated_selected.root");
  const TString perfAmbiFile = base + pathToFilesFromOutput
                         + (useIterativeTracking ?
                            "/performance_merged_ambi_tracks.root" :
                            "/performance_finding_ambi.root");
  const TString perfSeedFile = base + pathToFilesFromOutput
                         + (useIterativeTracking ?
                            "/performance_merged_seed.root" :
                            "/performance_seeding.root");
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
  bool perfMatched;

  uint32_t seedPerfParticle;
  uint32_t seedPerfEventId;
  bool seedPerfMatched;

  // set other branch addresses later to avoid loading vectors for this step
  tPerf->SetBranchAddress("event_nr", &perfEventId);
  tPerf->SetBranchAddress("particle_id_particle", &perfParticle);
  tPerf->SetBranchAddress("matched", &perfMatched);
  tSeedPerf->SetBranchAddress("event_nr", &seedPerfEventId);
  tSeedPerf->SetBranchAddress("particle_id_particle", &seedPerfParticle);
  tSeedPerf->SetBranchAddress("matched", &seedPerfMatched);

  // each key in outer map is event
  // inner map is (particle_id) ->  matched bool
  using MatchMap = std::unordered_map<unsigned, std::map<unsigned, bool>>;
  MatchMap matchMap(tPerf->GetEntries());
  MatchMap seedMatchMap(tSeedPerf->GetEntries());

  std::cout << "Loading matchingdetails from " << perfAmbiFile <<
            " with " << tPerf->GetEntries() << " entries." << std::endl;

  std::cout << "matchMap size: " << matchMap.size() << std::endl;

  for (Long64_t i = 0; i < tPerf->GetEntries(); i++) {
    tPerf->GetEntry(i);
    // check if event id already exists in map, if not create new entry
    if (matchMap.find(perfEventId) == matchMap.end()) {
      matchMap[perfEventId] = {};
    }
    matchMap[perfEventId][perfParticle] = perfMatched;
  }
  std::cout << "Loaded " << matchMap.size()
            << " entries from ambi matchingdetails" << std::endl;

  std::cout << "Loading matchingdetails from " << perfSeedFile <<
            " with " << tSeedPerf->GetEntries() << " entries." << std::endl;
  for (Long64_t i = 0; i < tSeedPerf->GetEntries(); i++) {
    tSeedPerf->GetEntry(i);
    if (seedMatchMap.find(seedPerfEventId) == seedMatchMap.end()) {
      seedMatchMap[seedPerfEventId] = {};
    }
    seedMatchMap[seedPerfEventId][seedPerfParticle] = seedPerfMatched;
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
  std::vector<bool> outMatched;
  std::vector<bool> outSeedMatched;

  tSim->SetBranchAddress("event_id", &simEventId);
  tSim->SetBranchAddress("particle", &simParticle);
  tSim->SetBranchAddress("generation", &simGeneration);
  tOut->Branch("matched", &outMatched);
  tOut->Branch("seed_matched", &outSeedMatched);

  Long64_t nSim   = tSim->GetEntries();
  Long64_t nFound = 0;
  Long64_t nSeedFound = 0;
  Long64_t nParticlesTotal = 0;

  for (Long64_t i = 0; i < nSim; i++) {
    tSim->GetEntry(i);
    if (i % 1000 == 0 && i) {
      std::cout << "Processed " << i << " / " << nSim
                << " events. Found matches for " << nFound
                << " particles." << std::endl;
    }
    if (simParticle->empty()) {  // There are no particles in this event
      tOut->Fill();  // fill the empty vectors, as we need one per event
      continue;
    }
    // set matched output vectors to default false for each particle
    outMatched.resize(simParticle->size(), false);
    outSeedMatched.resize(simParticle->size(), false);
    auto itEvent = matchMap.find(simEventId);

    // check each particle in event only if event has matching entry
    if (itEvent != matchMap.end()) {
      for (unsigned i_part = 0; i_part < simParticle->size(); i_part++) {
        // If only checking primaries, skip particle if generation is not 0
        if (checkPrimaryOnly && (*simGeneration)[i_part] != 0) {
          continue;
        }
        nParticlesTotal++;
        
        unsigned simPID = (*simParticle)[i_part];
        auto itPID = itEvent->second.find(simPID);
        // this simPID has an entry? sanity check since they all should have
        if (itPID != itEvent->second.end()) {
          outMatched[i_part] = itPID->second;
          nFound++;
        }
      }
    }

    auto itSeedEvent = seedMatchMap.find(simEventId);
    // event has matched seeds
    if (itSeedEvent != seedMatchMap.end()) {
      for (unsigned i_part = 0; i_part < simParticle->size(); i_part++) {
        if (checkPrimaryOnly && (*simGeneration)[i_part] != 0) {
          continue;
        }
        unsigned simPID = (*simParticle)[i_part];
        auto itSeedPID = itSeedEvent->second.find(simPID);
        // this simPID has an entry? sanity check since they all should have
        if (itSeedPID != itSeedEvent->second.end()) {
          outSeedMatched[i_part] = itSeedPID->second;
          nSeedFound++;
        }
      }
    }
    tOut->Fill();
    // empty the vectors for next event
    outMatched.clear();
    outSeedMatched.clear();
  }

  fOut->Write();

  fPerf->Close();
  fSeedPerf->Close();
  fOut->Close();
  fSim->Close();

  std::cout << "Matched " << nFound << " / " << nParticlesTotal << " particles." << std::endl;
  std::cout << "Seed matched " << nSeedFound << " / " << nParticlesTotal << " particles." << std::endl;
  std::cout << "Saved: " << outFile << std::endl;
}
