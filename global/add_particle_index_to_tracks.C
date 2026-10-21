#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "../Utils/BasicFunctions.h"

/*
 * Writes, per track, particle_index for the track summary: the global index of the
 * track's majority particle in particleFile
 * This index is -1 if the particle is not in particleFile. Filled for every track.
 * Output: friend file <trackFile without .root>_particle_index.root next to the input,
 * holding only the tree "trackindex" with the branch particle_index, one entry per
 * tracksummary entry (dCache files cannot be updated in place). Read it with
 * tracksummary->AddFriend("trackindex", <friend file>).
 */
void add_particle_index_to_tracks(const TString pathToFilesFromOutput,
                                  const bool onSTBC = false,
                                  const TString trackFile = "tracksummary_ambi-tracks-merged.root",
                                  const TString particleFile = "particles_digitized_selected.root") {
  // event, vertex_primary, vertex_secondary, particle, generation, sub_particle
  using Key = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;
  const TString dir = Utils::outputPathBase(onSTBC) + pathToFilesFromOutput + "/";
  const TString outFile = dir + TString(trackFile).ReplaceAll(".root", "_particle_index.root");

  TFile* fParticles = TFile::Open(dir + particleFile);
  TFile* fTracks = TFile::Open(dir + trackFile);
  TTree* tParticles = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
  TTree* tTracks = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
  if (!tParticles || !tTracks) {
    std::cerr << "Cannot read particles or tracksummary in " << dir << std::endl;
    return;
  }

  // global index of every particle
  unsigned event;
  std::vector<unsigned> *v1 = nullptr, *v2 = nullptr, *p = nullptr, *g = nullptr, *sp = nullptr;
  tParticles->SetBranchStatus("*", false);
  for (const char* name : {"event_id", "vertex_primary", "vertex_secondary", "particle", "generation", "sub_particle"})
    tParticles->SetBranchStatus(name, true);
  tParticles->SetBranchAddress("event_id", &event);
  tParticles->SetBranchAddress("vertex_primary", &v1);
  tParticles->SetBranchAddress("vertex_secondary", &v2);
  tParticles->SetBranchAddress("particle", &p);
  tParticles->SetBranchAddress("generation", &g);
  tParticles->SetBranchAddress("sub_particle", &sp);
  std::map<Key, int> indexOf;
  int nParticles = 0;
  for (Long64_t i = 0; i < tParticles->GetEntries(); i++) {
    tParticles->GetEntry(i);
    for (size_t k = 0; k < p->size(); k++)
      indexOf[{event, v1->at(k), v2->at(k), p->at(k), g->at(k), sp->at(k)}] = nParticles++;
  }

  // majority particle per track
  tTracks->SetBranchStatus("*", false);
  for (const char* name : {"event_nr", "majorityParticleId_vertex_primary", "majorityParticleId_vertex_secondary",
                           "majorityParticleId_particle", "majorityParticleId_generation",
                           "majorityParticleId_sub_particle"})
    tTracks->SetBranchStatus(name, true);
  tTracks->SetBranchAddress("event_nr", &event);
  tTracks->SetBranchAddress("majorityParticleId_vertex_primary", &v1);
  tTracks->SetBranchAddress("majorityParticleId_vertex_secondary", &v2);
  tTracks->SetBranchAddress("majorityParticleId_particle", &p);
  tTracks->SetBranchAddress("majorityParticleId_generation", &g);
  tTracks->SetBranchAddress("majorityParticleId_sub_particle", &sp);

  TFile* fOut = TFile::Open(outFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Cannot create " << outFile << std::endl;
    return;
  }
  TTree* tOut = new TTree("trackindex", "Global index of the majority particle in " + particleFile);
  std::vector<int> particleIndex;
  tOut->Branch("particle_index", &particleIndex);

  long nTracks = 0, nFound = 0;
  for (Long64_t i = 0; i < tTracks->GetEntries(); i++) {
    tTracks->GetEntry(i);
    particleIndex.assign(p->size(), -1);
    for (size_t k = 0; k < p->size(); k++) {
      auto it = indexOf.find({event, v1->at(k), v2->at(k), p->at(k), g->at(k), sp->at(k)});
      if (it != indexOf.end()) {
        particleIndex[k] = it->second;
        nFound++;
      }
    }
    nTracks += p->size();
    tOut->Fill();
  }
  tOut->Write();
  const bool failed = fOut->TestBit(TFile::kWriteError);
  fOut->Close();
  fTracks->Close();
  fParticles->Close();
  if (failed) {
    std::cerr << "Error writing " << outFile << ", remove it and rerun." << std::endl;
    return;
  }
  std::cout << "Majority particle of " << nFound << " / " << nTracks << " tracks found among the "
            << nParticles << " particles in " << particleFile << ". Saved: " << outFile << std::endl;
}
