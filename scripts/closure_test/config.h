#ifndef CLOSURE_TEST_CONFIG_H
#define CLOSURE_TEST_CONFIG_H

///////////////////////////////////////////////////////////////////////////////
// config.h
//
// Shared configuration for the MC23e cell-energy reweighting closure test.
// Defines branch names, binning, geometry, distortion parameters, and helper
// functions adapted for MC23e ntuples.
//
// Based on supervisor's var.h and tree.h from showershapereweighting package.
// Reference: ATL-COM-PHYS-2021-640, Section 5.2 (cell reweighting method)
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

namespace closure {

    // ======================================================================
    // Cluster geometry (same as supervisor's var.h)
    // ======================================================================
    const int kEtaSize     = 7;
    const int kPhiSize     = 11;
    const int kClusterSize = kEtaSize * kPhiSize;   // 77

    // Central cell index (eta=3, phi=5) → flat index = 3*11+5 = 38
    const int kCentralCell = 3 * kPhiSize + 5;       // 38

    // Layer 2 cell sizes
    const double kDeltaEta = 0.025;
    const double kDeltaPhi = 0.0245;

    // ======================================================================
    // Eta binning (14 bins, same as supervisor's var.h)
    // ======================================================================
    const int kNEtaBins = 14;
    const double kEtaLimits[kNEtaBins + 1] = {
        0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.37,
        1.52, 1.6, 1.80, 2.0, 2.2, 2.4
    };

    // Energy fraction bins for CellCorrection TProfile
    // (matching supervisor's EBins in tree.h — NOT var.h which is for ClusterCorrection in GeV)
    // Cell fractions f_k = e_k / E_tot live in [~-0.1, ~0.5], so these bins
    // give fine resolution where data actually populates.
    const int kNFracBins = 12;
    const double kFracBins[kNFracBins + 1] = {
        -1., -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 1.
    };

    // ======================================================================
    // Branch names for MC23e ntuples
    // ======================================================================
    const char* const kTreeName       = "tree";
    const char* const kCellBranch     = "photon.7x11ClusterLr2E";
    const char* const kCellSizeBranch = "photon.7x11ClusterLr2Size";
    const char* const kEtaBranch      = "photon_cluster.eta2";

    // ======================================================================
    // Systematic distortion pattern for pseudo-data generation
    // ======================================================================
    // Mimics typical data-MC differences: data has narrower showers, so
    // central cells are enhanced and outer cells suppressed.
    // Returns a base bias in [-1, +1] that is scaled by distortionLevel.
    //
    // Pattern:  center -> +1, ring1 -> +0.4, ring2 -> -0.4, outer -> -1.0
    // Usage:    e'_k = e_k * (1 + level*getBiasPattern(k) + Gaus(0, noiseSigma))
    // ======================================================================
    inline double getBiasPattern(int cellIdx) {
        int eta_row = cellIdx / kPhiSize;     // 0-6
        int phi_col = cellIdx % kPhiSize;     // 0-10
        int deta = std::abs(eta_row - 3);
        int dphi = std::abs(phi_col - 5);
        int dist = std::max(deta, dphi);      // Chebyshev distance from center

        if (dist == 0) return +1.0;   // central cell
        if (dist == 1) return +0.4;   // first ring
        if (dist == 2) return -0.4;   // second ring
        return -1.0;                  // outer cells
    }

    // ======================================================================
    // Cluster quality checks for MC23e ntuples
    // ======================================================================
    // The supervisor's ntuples are pre-selected matched data-MC pairs with
    // tight cuts. MC23e ntuples contain all photon candidates including
    // low-pT, miscentered, and noise-dominated clusters. These functions
    // replicate the implicit quality cuts of the supervisor's selection.

    // Clamp negative cell energies to zero.
    // Negative energies come from electronic noise subtraction and create
    // unphysical fractions that corrupt TProfile corrections.
    inline void clampCellEnergies(std::vector<double>& cells) {
        for (int k = 0; k < kClusterSize; ++k) {
            if (cells[k] < 0) cells[k] = 0;
        }
    }

    // Check that the central cell has the highest energy fraction.
    // Miscentered clusters (where the hottest cell is off-center) produce
    // anomalous fractions that bias the correction derivation.
    inline bool isHealthyCluster(const std::vector<double>& cells) {
        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k)
            Etot += cells[k];
        if (Etot <= 0) return false;

        double centralFrac = cells[kCentralCell] / Etot;
        for (int k = 0; k < kClusterSize; ++k) {
            if (k == kCentralCell) continue;
            if (cells[k] / Etot > centralFrac) return false;
        }
        return true;
    }

    // ======================================================================
    // Find eta bin index for a given |eta| value
    // Returns -1 if outside range
    // ======================================================================
    inline int findEtaBin(double abseta) {
        for (int n = 0; n < kNEtaBins; ++n) {
            if (abseta >= kEtaLimits[n] && abseta < kEtaLimits[n + 1])
                return n;
        }
        return -1;
    }

    // ======================================================================
    // Energy in a sub-window of the 7x11 cluster
    // (same logic as supervisor's Energy() in var.h)
    //
    // eta, phi: window dimensions (must be odd)
    // Returns sum of energies in the centred window
    // ======================================================================
    inline double Energy(int eta, int phi,
                         const std::vector<double>& clusterEnergy) {
        int etaMin = kEtaSize - (kEtaSize + eta) / 2;
        int etaMax = kEtaSize - (kEtaSize - eta) / 2;
        int phiMin = kPhiSize - (kPhiSize + phi) / 2;
        int phiMax = kPhiSize - (kPhiSize - phi) / 2;
        double sumE = 0;
        for (int e = etaMin; e < etaMax; ++e) {
            for (int p = phiMin; p < phiMax; ++p) {
                sumE += clusterEnergy[p + kPhiSize * e];
            }
        }
        return sumE;
    }

    // ======================================================================
    // Compute R_eta = E(3x7) / E(7x7)
    // ======================================================================
    inline double calcReta(const std::vector<double>& cells) {
        double E7x7 = Energy(7, 7, cells);
        double E3x7 = Energy(3, 7, cells);
        if (E7x7 <= 0) return -999.0;
        return E3x7 / E7x7;
    }

    // ======================================================================
    // Compute R_phi = E(3x3) / E(3x7)
    // ======================================================================
    inline double calcRphi(const std::vector<double>& cells) {
        double E3x7 = Energy(3, 7, cells);
        double E3x3 = Energy(3, 3, cells);
        if (E3x7 <= 0) return -999.0;
        return E3x3 / E3x7;
    }

    // ======================================================================
    // Compute w_eta_2 from cell energies using relative eta positions
    //
    // w_eta_2 = sqrt( sum(E_i * eta_i^2)/sum(E_i)
    //               - (sum(E_i * eta_i)/sum(E_i))^2 )
    //
    // Computed in a 3(eta) x 5(phi) window centred on the cluster.
    // Uses relative eta positions: (-1, 0, +1) * kDeltaEta
    //
    // NOTE: Supervisor's calcWeta() in var.h uses absolute eta from
    // clusterEta branch. MC23e ntuples lack this branch, so we use
    // relative positions instead. Results are equivalent up to a
    // constant offset that cancels in the width calculation.
    // ======================================================================
    inline double calcWeta2(const std::vector<double>& cells) {
        const int neta_w = 3, nphi_w = 5;
        int eta_start = kEtaSize / 2 - neta_w / 2;   // 2
        int phi_start = kPhiSize / 2 - nphi_w / 2;    // 3

        double sumE = 0, sumEeta = 0, sumEeta2 = 0;
        for (int ie = 0; ie < neta_w; ++ie) {
            double eta_rel = (ie - neta_w / 2) * kDeltaEta;
            for (int ip = 0; ip < nphi_w; ++ip) {
                int idx = (eta_start + ie) * kPhiSize + (phi_start + ip);
                double e = cells[idx];
                sumE     += e;
                sumEeta  += e * eta_rel;
                sumEeta2 += e * eta_rel * eta_rel;
            }
        }
        if (sumE <= 0) return -999.0;
        double var = sumEeta2 / sumE - std::pow(sumEeta / sumE, 2);
        if (var < 0) return -999.0;
        return std::sqrt(var);
    }

}  // namespace closure

#endif // CLOSURE_TEST_CONFIG_H
