// Acceptance.h
//
// This file contains utility functions and classes related to particle acceptance criteria.
#pragma once

#include <vector>
#include <unordered_set>
#include <iostream>
#include <limits>
#include <algorithm>

#include <TString.h>

namespace Utils {
  class Acceptance {
  public:
    // only create instance and then add ranges and particles to accept
    Acceptance()
      : m_pt_range(nullptr),
        m_eta_range(nullptr),
        m_abs_eta_range(nullptr),
        m_phi_range(nullptr) {}

    ~Acceptance() = default;

    // make add functions chainable
    Acceptance& addPtRange(float ptmin, float ptmax) {
      m_pt_range = new std::tuple<float, float>(ptmin, ptmax);
      return *this;
    }

    Acceptance& addEtaRange(float etamin, float etamax) {
      if (m_abs_eta_range) {
        std::cout << "Cannot add eta range: absolute eta range is already set."
                  << std::endl;
      } else {
        m_eta_range = new std::tuple<float, float>(etamin, etamax);
      }
      return *this;
    }

    Acceptance& addAbsEtaRange(float abs_etamin, float abs_etamax) {
      if (m_eta_range) {
        std::cout << "Cannot add absolute eta range: eta range is already set."
                  << std::endl;
      } else if (abs_etamin < 0 || abs_etamax < 0) {
        std::cout << "Cannot add absolute eta range: values must be non-negative."
                  << std::endl;
      } else {
        m_abs_eta_range = new std::tuple<float, float>(abs_etamin, abs_etamax);
      }
      return *this;
    }

    Acceptance& addPhiRange(float phimin, float phimax) {
      m_phi_range = new std::tuple<float, float>(phimin, phimax);
      return *this;
    }

    Acceptance& addParticleToAccept(int particle) {
      // check if particles to reject exists
      if (!m_particles_to_reject.empty()) {
        std::cout << "Cannot add particle to accept: particles work on per-reject basis."
                  << std::endl;
      } else {
        m_particles_to_accept.insert(particle);
      }
      return *this;
    }

    Acceptance& addParticleToReject(int particle) {
      // check if particles to accept exists
      if (!m_particles_to_accept.empty()) {
        std::cout << "Cannot add particle to reject: particles work on per-accept basis."
                  << std::endl;
      } else {
        m_particles_to_reject.insert(particle);
      }
      return *this;
    }

    Acceptance& clearParticlesToAccept() {
      m_particles_to_accept.clear();
      return *this;
    }

    Acceptance& clearParticlesToReject() {
      m_particles_to_reject.clear();
      return *this;
    }

    Acceptance& standardAcceptance(const TString& type="ALICE3") {
      if (type == "ALICE3") {
        addPtRange(0.5, 10.0);
        addEtaRange(-2.5, 2.5);
      } else if (type == "loopers_central") {
        addPtRange(0.01, 0.2);
        addEtaRange(-1.0, 1.0);
      } else if (type == "OT_discs_coverage") {
        addPtRange(0.5, 10.0);
        addAbsEtaRange(1.5, 2.5);
      } else {
        std::cout << "Unknown standard acceptance type: " << type << std::endl;
      }
      return *this;
    }

    bool isAccepted(float pt, float eta, float phi, int particle) const {
      // first check if range is set, then whether particle falls within range
      if (m_pt_range && (pt < std::get<0>(*m_pt_range) || pt > std::get<1>(*m_pt_range))) {
        return false;
      }
      if (m_eta_range && (eta < std::get<0>(*m_eta_range) || eta > std::get<1>(*m_eta_range))) {
        return false;
      }
      if (m_abs_eta_range && (std::abs(eta) < std::get<0>(*m_abs_eta_range) ||
                              std::abs(eta) > std::get<1>(*m_abs_eta_range))) {
        return false;
      }
      if (m_phi_range && (phi < std::get<0>(*m_phi_range) || phi > std::get<1>(*m_phi_range))) {
        return false;
      }
      if (!m_particles_to_accept.empty() && m_particles_to_accept.find(particle) == m_particles_to_accept.end()) {
        return false;
      }
      if (!m_particles_to_reject.empty() && m_particles_to_reject.find(particle) != m_particles_to_reject.end()) {
        return false;
      }
      return true;
    }

  private:
    // pT in GeV
    std::tuple<float, float>* m_pt_range;
    std::tuple<float, float>* m_eta_range;
    std::tuple<float, float>* m_abs_eta_range;
    std::tuple<float, float>* m_phi_range;

    std::unordered_set<int> m_particles_to_accept;
    std::unordered_set<int> m_particles_to_reject;
  };
}  // namespace Utils