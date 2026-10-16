#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "../Utils/BasicFunctions.h"

/*
 * Adds, per track, particle_index to the track summary: the global index of the
 * track's majority particle in particleFile
 * This index is -1 if the particle is not in particleFile. Filled for every track.
 * Output: Simply adds to the input track file.
 */
void add_particle_index_to_tracks(const TString pathToFilesFromOutput,
                                  const bool onSTBC = false,
                                  const TString trackFile = "tracksummary_ambi-tracks-merged.root",
                                  const TString particleFile = "particles_digitized_selected.root") {
  // event, vertex_primary, vertex_secondary, particle, generation, sub_particle
  using Key = std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned, unsigned>;
  const TString dir = Utils::outputPathBase(onSTBC) + pathToFilesFromOutput + "/";

  TFile* fParticles = TFile::Open(dir + particleFile);
  TFile* fTracks = TFile::Open(dir + trackFile, "UPDATE");
  TTree* tParticles = fParticles && !fParticles->IsZombie() ? fParticles->Get<TTree>("particles") : nullptr;
  TTree* tTracks = fTracks && !fTracks->IsZombie() ? fTracks->Get<TTree>("tracksummary") : nullptr;
  if (!tParticles || !tTracks) {
    std::cerr << "Cannot read particles or tracksummary in " << dir << std::endl;
    return;
  }
  if (tTracks->GetBranch("particle_index")) {
    std::cout << "particle_index already in " << dir + trackFile << ", nothing to do." << std::endl;
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
  tTracks->SetBranchAddress("event_nr", &event);
  tTracks->SetBranchAddress("majorityParticleId_vertex_primary", &v1);
  tTracks->SetBranchAddress("majorityParticleId_vertex_secondary", &v2);
  tTracks->SetBranchAddress("majorityParticleId_particle", &p);
  tTracks->SetBranchAddress("majorityParticleId_generation", &g);
  tTracks->SetBranchAddress("majorityParticleId_sub_particle", &sp);

  std::vector<int> particleIndex;
  TBranch* bIndex = tTracks->Branch("particle_index", &particleIndex);

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
    bIndex->Fill();
  }
  fTracks->cd();
  tTracks->Write("", TObject::kOverwrite);
  fTracks->Close();
  fParticles->Close();
  std::cout << "Majority particle of " << nFound << " / " << nTracks << " tracks found among the "
            << nParticles << " particles in " << particleFile << ". Added to " << dir + trackFile << std::endl;
}
