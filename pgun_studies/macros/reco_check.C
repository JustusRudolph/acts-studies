// This is a super simple macro to check which particles were reconstructed
#include <TFile.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <iostream>
#include <vector>

void print_particle_info_from_event(unsigned i_event, TDatabasePDG* pdg, 
                                    std::vector<bool>* matched,
                                    std::vector<int>* pid,
                                    std::vector<unsigned>* particle_barcode,
                                    std::vector<unsigned>* subparticle_barcode,
                                    std::vector<unsigned>* generation,
                                    std::vector<unsigned>* nHits,
                                    std::vector<float>* p,
                                    std::vector<float>* vx,
                                    std::vector<float>* vy,
                                    std::vector<float>* vz,
                                    std::vector<float>* vt,
                                    std::vector<float>* px,
                                    std::vector<float>* py,
                                    std::vector<float>* pz) {

  for (size_t i_part = 0; i_part < particle_barcode->size(); i_part++) {
    // get radial component of production vertex of each part/subpart
    float vr = std::sqrt(vx->at(i_part)*vx->at(i_part)
                        + vy->at(i_part)*vy->at(i_part));
    // \dot{r} = (\vec{r} . \vec {v}) / |r|
    float pr = (px->at(i_part)*vx->at(i_part)
              + py->at(i_part)*vy->at(i_part)) / vr;
    std::cout << "\tParticle " << particle_barcode->at(i_part)
              << " has subparticle " << subparticle_barcode->at(i_part)
              << " of generation " << generation->at(i_part)
              << " which is a " << pdg->GetParticle(pid->at(i_part))->GetName()
              << " with momentum " << p->at(i_part) << " and " << nHits->at(i_part)
              << " hits. Production vertex (r, z, t) = (" << vr << ", "
              << vz->at(i_part) << ", " << vt->at(i_part) << ") in direction ("
              << pr << ", "  << pz->at(i_part) << ")" << std::endl;
  }
}

void print_if_subparticle(unsigned i_event, TDatabasePDG* pdg, 
                          std::vector<bool>* matched,
                          std::vector<int>* pid,
                          std::vector<unsigned>* particle_barcode,
                          std::vector<unsigned>* subparticle_barcode,
                          std::vector<unsigned>* generation,
                          std::vector<unsigned>* nHits,
                          std::vector<float>* p,
                          std::vector<float>* vx,
                          std::vector<float>* vy,
                          std::vector<float>* vz,
                          std::vector<float>* vt,
                          std::vector<float>* px,
                          std::vector<float>* py,
                          std::vector<float>* pz) {
  // check if there is a subparticle (so >1)
  if (std::find_if(subparticle_barcode->begin(), subparticle_barcode->end(),
            [](unsigned code) { return code > 1; }) != subparticle_barcode->end()) {
    
    std::cout << "In event " << i_event << " there is a decay chain ";
    if (matched->size()) {
      std::cout << " whose primary particle is " << (matched->front() ? "" : "NOT ")
                << "matched." << std::endl;
    } else {
      std::cout << " but no matching details are available." << std::endl;
    }
    print_particle_info_from_event(i_event, pdg, matched, pid,
                                   particle_barcode, subparticle_barcode,
                                   generation, nHits, p, vx, vy, vz, vt,
                                   px, py, pz);
  }
}

void print_if_no_matching_details(unsigned i_event, TDatabasePDG* pdg, 
                                   std::vector<bool>* matched,
                                   std::vector<int>* pid,
                                   std::vector<unsigned>* particle_barcode,
                                   std::vector<unsigned>* subparticle_barcode,
                                   std::vector<unsigned>* generation,
                                   std::vector<unsigned>* nHits,
                                   std::vector<float>* p,
                                   std::vector<float>* vx,
                                   std::vector<float>* vy,
                                   std::vector<float>* vz,
                                   std::vector<float>* vt,
                                   std::vector<float>* px,
                                   std::vector<float>* py,
                                   std::vector<float>* pz) {
  if (matched->empty()) {
    std::cout << "In event " << i_event << ", no matching details are available." << std::endl;
    print_particle_info_from_event(i_event, pdg, matched, pid,
                                   particle_barcode, subparticle_barcode,
                                   generation, nHits, p, vx, vy, vz, vt,
                                   px, py, pz);
  }
}

void reco_check(TString input_base_path = "../output_tutorial/test_electron/") {
  TDatabasePDG* pdg = TDatabasePDG::Instance();

  TFile* perf_file = TFile::Open(input_base_path + "performance_finding_ambi.root");
  TFile* sim_file = TFile::Open(input_base_path + "particles_simulation.root");
    
  if (!perf_file || perf_file->IsZombie()) {
    std::cout << "Warning: Could not open " << input_base_path + "performance_finding_ambi.root" << std::endl;
    return;
  }
  if (!sim_file || sim_file->IsZombie()) {
    std::cout << "Warning: Could not open " << input_base_path + "particles_simulation.root" << std::endl;
    return;
  }

  TTree* matchingdetails = (TTree*) perf_file->Get("matchingdetails");
  TTree* particles_tree = (TTree*) sim_file->Get("particles");
  
  std::vector<bool>* matched = nullptr;
  matchingdetails->SetBranchAddress("matched", &matched);

  std::cout << "Number of events in matchingdetails tree: " << matchingdetails->GetEntries() << std::endl;
  for (unsigned i_ev = 0; i_ev < matchingdetails->GetEntries(); i_ev++) {
      matchingdetails->GetEntry(i_ev);
      for (size_t i_mcp = 0; i_mcp < matched->size(); i_mcp++) {
          if (matched->at(i_mcp)) {
              std::cout << "Event " << i_ev << " MCP " << i_mcp << " was matched!" << std::endl;
          } else {
              std::cout << "Event " << i_ev << " MCP " << i_mcp << " was NOT matched!" << std::endl;
          }
      }
  }

  // check what the subparticles are
  std::vector<int>* pid{nullptr};
  std::vector<unsigned>* particle_barcode{nullptr}, *subparticle_barcode{nullptr},
                        *nHits{nullptr}, *generation{nullptr};
  std::vector<float>* p{nullptr}, *vx{nullptr}, *vy{nullptr}, *vz{nullptr},
                      *vt{nullptr}, *px{nullptr}, *py{nullptr}, *pz{nullptr};

  particles_tree->SetBranchAddress("particle_type", &pid);
  particles_tree->SetBranchAddress("particle", &particle_barcode);
  particles_tree->SetBranchAddress("sub_particle", &subparticle_barcode);
  particles_tree->SetBranchAddress("generation", &generation);
  particles_tree->SetBranchAddress("number_of_hits", &nHits);
  particles_tree->SetBranchAddress("p", &p);
  particles_tree->SetBranchAddress("vx", &vx);
  particles_tree->SetBranchAddress("vy", &vy);
  particles_tree->SetBranchAddress("vz", &vz);
  particles_tree->SetBranchAddress("vt", &vt);
  particles_tree->SetBranchAddress("px", &px);
  particles_tree->SetBranchAddress("py", &py);
  particles_tree->SetBranchAddress("pz", &pz);

  // create histogram of number of hits of primary particles
  TH1F* h_nHits_primary = new TH1F("h_nHits_primary", "Number of Hits of Primary Particles", 15, 0, 15);

  std::cout << "Number of events in particles tree: " << particles_tree->GetEntries() << std::endl;
  for (unsigned i_event = 0; i_event < particles_tree->GetEntries(); i_event++) {
    particles_tree->GetEntry(i_event);
    matchingdetails->GetEntry(i_event);

    if (particle_barcode->size() != subparticle_barcode->size()) {
      std::cout << "In event " << i_event << ", particle and subparticle sizes differ: "
                << particle_barcode->size() << " vs " << subparticle_barcode->size()
                << std::endl;
      continue;
    }
    h_nHits_primary->Fill(nHits->front());  // primary particle is first entry

    print_if_no_matching_details(i_event, pdg, matched, pid,
                                 particle_barcode, subparticle_barcode,
                                 generation, nHits, p, vx, vy, vz, vt,
                                 px, py, pz);
  }
  // write out hist
  TFile* output_file = TFile::Open("hists/hits/reco_check.root", "RECREATE");
  h_nHits_primary->Write();
  output_file->Close();
}