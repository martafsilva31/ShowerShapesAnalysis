#ifndef DATA_MC_CONFIG_H
#define DATA_MC_CONFIG_H

///////////////////////////////////////////////////////////////////////////////
// config.h
//
// Shared configuration for the data-MC comparison and cell-energy reweighting
// pipeline.  Extends the closure_test/config.h conventions with offline
// selection cuts from ATL-COM-PHYS-2025-662 (Section 3).
//
// Data files:  egam3.root  (Z→eeγ data),    egam4.root  (Z→μμγ data)
// MC files:    mc_eegamma.root,              mc_mumugamma.root
//
// Ntuple-level cuts already applied by NTupleMaker (Zeeg/Zmumug):
//   - Trigger, GRL, detector quality, primary vertex
//   - ≥1 photon + ≥2 leptons with PassSelection
//   - Opposite-charge, trigger-matched leptons
//   - ΔR(γ,e) > 0.2 / ΔR(γ,μ) > 0.01
//   - 40 < mll < 120 GeV,  45 < mllg < 125 GeV
//
// Additional offline cuts (applied here):
//   - FSR mass windows:  40 < mll < 83 GeV,  80 < mllg < 100 GeV
//   - Photon |η| < 2.37, excluding crack 1.37–1.52
//   - Photon tight ID (tightOffInc)       [toggle]
//   - FixedCutTight isolation (isotight)   [toggle]
//   - Cell quality: 77 cells, central cell hottest  [for cell-level studies]
//   - ΔR(γ,l) > 0.4                                [for cell reweighting]
//   - Truth match (MC only)                         [toggle]
//
// Reference: ATL-COM-PHYS-2025-662, ATL-COM-PHYS-2021-640
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include <TChain.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

namespace datamc {

    // ======================================================================
    // Cluster geometry (same as closure_test/config.h)
    // ======================================================================
    const int kEtaSize     = 7;
    const int kPhiSize     = 11;
    const int kClusterSize = kEtaSize * kPhiSize;   // 77
    const int kCentralCell = 3 * kPhiSize + 5;       // 38
    const double kDeltaEta = 0.025;
    const double kDeltaPhi = 0.0245;

    // ======================================================================
    // FSR selection mass windows (from note Section 3.3)
    // The NTupleMaker uses wider windows (mll: 40-120, mllg: 45-125);
    // we tighten them here to isolate radiative Z decays.
    // ======================================================================
    const double kMllMin  = 40.0;   // GeV
    const double kMllMax  = 83.0;   // GeV
    const double kMllgMin = 80.0;   // GeV
    const double kMllgMax = 100.0;  // GeV

    // ======================================================================
    // Photon pT minimum (note Section 3.1)
    // ======================================================================
    const double kPhotonPtMin = 7.0;    // GeV

    // ======================================================================
    // Photon eta acceptance (note Section 3.1)
    // ======================================================================
    const double kEtaMax       = 2.37;
    const double kCrackEtaMin  = 1.37;
    const double kCrackEtaMax  = 1.52;

    // ======================================================================
    // Overlap removal for cell reweighting (note Section 3.5)
    // ======================================================================
    const double kDRPhotonLeptonMin = 0.4;

    // Lepton-lepton separation (note Section 3.2)
    const double kDRLeptonLeptonMin = 0.2;

    // ======================================================================
    // Eta binning (14 bins, same as supervisor's var.h and closure_test)
    // ======================================================================
    const int kNEtaBins = 14;
    const double kEtaLimits[kNEtaBins + 1] = {
        0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.37,
        1.52, 1.6, 1.80, 2.0, 2.2, 2.4
    };

    // ======================================================================
    // pT binning for corrections (note Section 5.3)
    // ======================================================================
    const int kNPtBins = 6;
    const double kPtLimits[kNPtBins + 1] = {
        10, 15, 20, 25, 30, 40, 1000
    };

    // Energy fraction bins for TProfile (Francisco's tree.h non-uniform bins)
    // Full range [-1, 1] to capture all cell fractions including
    // negative values from noise cells.
    const int kNFracBins = 12;
    const double kFracBins[kNFracBins + 1] = {
        -1., -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 1.
    };

    // Sigma ratio cap for M3 shift+stretch — prevents extreme corrections
    // in low-statistics eta bins (e.g. bins 5 and 9 had ratios up to 67×).
    const double kSigmaRatioMin = 0.5;
    const double kSigmaRatioMax = 2.0;

    // ======================================================================
    // Branch names for MC23e / data24 ntuples
    // ======================================================================
    const char* const kTreeName       = "tree";
    const char* const kCellBranch     = "photon.7x11ClusterLr2E";
    const char* const kCellEtaBranch  = "photon.7x11ClusterLr2Eta";
    const char* const kCellPhiBranch  = "photon.7x11ClusterLr2Phi";
    const char* const kCellSizeBranch = "photon.7x11ClusterLr2Size";
    const char* const kEtaBranch      = "photon_cluster.eta2";
    const char* const kPhotonPtBranch = "photon.pt";
    const char* const kPhotonEtaBranch = "photon.eta";
    const char* const kPhotonPhiBranch = "photon.phi";
    const char* const kPhotonEBranch   = "photon_cluster.e";
    const char* const kMllBranch      = "mll";
    const char* const kMllgBranch     = "mllg";
    const char* const kIsConvBranch   = "photon.isconv";
    const char* const kEtamax2Branch   = "photon_cluster.etamax2";

    // Photon ID and isolation
    const char* const kTightIDBranch  = "photon.tightOffInc";
    const char* const kIsoTightBranch = "photon.isotight";

    // Overlap removal
    const char* const kDRLepton1Branch = "lepton1.dr_ph";
    const char* const kDRLepton2Branch = "lepton2.dr_ph";

    // MC-specific
    const char* const kTruthMatchBranch = "photon.istruthmatch";
    const char* const kMCWeightBranch   = "event.mcwgt";
    const char* const kPUWeightBranch   = "event.muwgt";
    const char* const kBSWeightBranch   = "event.beamspotwgt";
    const char* const kLepton1SFBranch  = "lepton1.SF";
    const char* const kLepton2SFBranch  = "lepton2.SF";

    // Lepton separation
    const char* const kDRllBranch       = "dRll";

    // Unfudged shower shapes (raw from calorimeter)
    const char* const kRetaUnfBranch  = "photon.unfudged_reta";
    const char* const kRphiUnfBranch  = "photon.unfudged_rphi";
    const char* const kWeta2UnfBranch = "photon.unfudged_weta2";

    // Fudged shower shapes (after ATLAS corrections)
    const char* const kRetaFudBranch  = "photon.reta";
    const char* const kRphiFudBranch  = "photon.rphi";
    const char* const kWeta2FudBranch = "photon.weta2";

    // ======================================================================
    // Inline helpers
    // ======================================================================

    inline bool inCrack(double abseta) {
        return (abseta > kCrackEtaMin && abseta < kCrackEtaMax);
    }

    inline bool passEtaAcceptance(double abseta) {
        return (abseta < kEtaMax && !inCrack(abseta));
    }

    inline bool passFSRMassWindow(double mll, double mllg) {
        return (mll > kMllMin && mll < kMllMax &&
                mllg > kMllgMin && mllg < kMllgMax);
    }

    inline int findEtaBin(double abseta) {
        for (int n = 0; n < kNEtaBins; ++n) {
            if (abseta >= kEtaLimits[n] && abseta < kEtaLimits[n + 1])
                return n;
        }
        return -1;
    }

    inline int findPtBin(double pt) {
        for (int n = 0; n < kNPtBins; ++n) {
            if (pt >= kPtLimits[n] && pt < kPtLimits[n + 1])
                return n;
        }
        return -1;
    }

    // ======================================================================
    // Combined event selection for cell reweighting (note §3.5)
    // All values in GeV (NTupleMaker converts from MeV at write time).
    // ======================================================================
    // ======================================================================
    // Helper: create a TChain from a possibly comma-separated file list.
    // Single files work transparently (one Add call).
    // ======================================================================
    inline TChain* makeChain(const char* files) {
        TChain* chain = new TChain(kTreeName);
        TString flist(files);
        TObjArray* arr = flist.Tokenize(",");
        for (int i = 0; i < arr->GetEntries(); ++i) {
            TString path = ((TObjString*)arr->At(i))->GetString();
            path = path.Strip(TString::kBoth);
            chain->Add(path);
            std::cout << "  Chain: " << path << std::endl;
        }
        delete arr;
        return chain;
    }

    inline bool passSelection(double photonPt,
                              double clusterEta,
                              double mll,
                              double mllg,
                              double dR_l1_ph,
                              double dR_l2_ph,
                              double dR_ll,
                              int cellSize,
                              bool isTruthMatch,
                              bool isMC,
                              int isTight,
                              bool isConv) {
        // Tight photon ID (matching Francisco's Run 2 selection)
        if (isTight != 1) return false;
        // Unconverted photons only
        if (isConv) return false;
        if (photonPt < kPhotonPtMin) return false;
        double abseta = std::fabs(clusterEta);
        if (!passEtaAcceptance(abseta)) return false;
        if (!passFSRMassWindow(mll, mllg)) return false;
        if (dR_l1_ph < kDRPhotonLeptonMin) return false;
        if (dR_l2_ph < kDRPhotonLeptonMin) return false;
        if (dR_ll < kDRLeptonLeptonMin) return false;
        if (cellSize != kClusterSize) return false;
        if (isMC && !isTruthMatch) return false;
        return true;
    }

    // ======================================================================
    // Cluster quality
    //
    // isHealthyCluster: checks that central cell (k=38) carries the most
    // energy among all 77 cells, and total energy > 0.
    // Works on RAW cells (no clamping needed).  Matches Francisco's
    // NTUP.cxx quality cut: maxK == 38.
    // ======================================================================
    inline bool isHealthyCluster(const std::vector<double>& cells) {
        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k)
            Etot += cells[k];
        if (Etot <= 0) return false;

        // Central cell must be the hottest (by raw energy, not fraction)
        double centralE = cells[kCentralCell];
        for (int k = 0; k < kClusterSize; ++k) {
            if (k == kCentralCell) continue;
            if (cells[k] > centralE) return false;
        }
        return true;
    }

    // Legacy: clamp negative cell energies to zero.
    // NOT used in the reweighting pipeline (Francisco doesn't clamp).
    // Kept for backward compat with closure_test.
    inline void clampCellEnergies(std::vector<double>& cells) {
        for (int k = 0; k < kClusterSize; ++k) {
            if (cells[k] < 0) cells[k] = 0;
        }
    }

    // ======================================================================
    // Shower shape calculations — index-based (legacy)
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

    inline double calcReta(const std::vector<double>& cells) {
        double E7x7 = Energy(7, 7, cells);
        double E3x7 = Energy(3, 7, cells);
        if (E7x7 <= 0) return -999.0;
        return E3x7 / E7x7;
    }

    inline double calcRphi(const std::vector<double>& cells) {
        double E3x7 = Energy(3, 7, cells);
        double E3x3 = Energy(3, 3, cells);
        if (E3x7 <= 0) return -999.0;
        return E3x3 / E3x7;
    }

    inline double calcWeta2Raw(const std::vector<double>& cells) {
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

    // ======================================================================
    // ATLAS weta2 position-dependent correction (egammaqweta2c)
    //
    // The raw weta2 (std dev of eta-energy distribution in 3×5 window)
    // depends on where the shower lands within the calorimeter cell.
    // ATLAS subtracts a polynomial correction P(etarel) to remove this
    // dependence.  Coefficients from Athena egammaqweta2c.cxx.
    //
    // etarel = fmod(|eta - etacell - 0.025/2|, 0.025) / 0.025
    //   where eta    = cluster eta in sampling 2  (photon_cluster.eta2)
    //         etacell = eta of hottest cell in L2  (photon_cluster.etamax2)
    // ======================================================================
    inline double weta2RelPosition(float eta, float etacell) {
        double x = std::fabs(eta - etacell - 0.025 / 2.);
        return std::fmod(x, 0.025) / 0.025;
    }

    inline float weta2Correct(float eta, float etacell, float weta2_raw) {
        if (eta == -999.f || etacell == -999.f) return weta2_raw;
        float aeta = std::fabs(eta);
        float u = weta2RelPosition(eta, etacell);

        // Polynomial coefficients: correction = P0 + P1*u + P2*u^2
        static const float P0A[3] = { 0.0045f,   0.005375f, -0.0562f  };
        static const float P1A[3] = {-0.0016f,  -0.0215f,    0.114f   };
        static const float P2A[3] = {-0.0866f,   0.0215f,   -0.053f   };
        static const float P0B[3] = { 0.0039f,   0.005075f, -0.0324f  };
        static const float P1B[3] = { 0.00816f, -0.0203f,    0.0653f  };
        static const float P2B[3] = {-0.145f,    0.0203f,   -0.0286f  };
        static const float P0C[2] = { 0.0047f,   0.0035f  };
        static const float P1C[2] = {-0.0184f,  -0.0139f  };
        static const float P2C[2] = { 0.0180f,   0.0137f  };

        auto poly = [&](const float* P0, const float* P1, const float* P2,
                        int region) -> float {
            return weta2_raw - (P0[region] + P1[region]*u + P2[region]*u*u);
        };

        if (aeta < 0.8f) {
            int r = (u < 0.1f) ? 0 : (u < 0.9f) ? 1 : 2;
            return poly(P0A, P1A, P2A, r);
        }
        if (aeta < 1.8f) {
            int r = (u < 0.1f) ? 0 : (u < 0.9f) ? 1 : 2;
            return poly(P0B, P1B, P2B, r);
        }
        if (aeta < 2.0f) return poly(P0C, P1C, P2C, 0);
        if (aeta < 2.5f) return poly(P0C, P1C, P2C, 1);
        return weta2_raw;
    }

    // Corrected weta2: raw cell computation + ATLAS position correction.
    // This should match the branch value (photon.unfudged_weta2).
    inline double calcWeta2(const std::vector<double>& cells,
                            float eta, float etacell) {
        double raw = calcWeta2Raw(cells);
        if (raw < 0) return raw;
        return weta2Correct(eta, etacell, raw);
    }

    // ======================================================================
    // Shower shape calculations — position-based
    //
    // Use stored cell eta/phi positions to define geometric windows,
    // matching Athena's CaloCellList approach.  This avoids the grid
    // ordering issue at eta ≈ 0 (where CaloFillRectangularCluster's
    // output wraps across the eta = 0 boundary) and the non-uniform
    // cell geometry at the barrel-endcap transition (|eta| ≈ 1.3-1.5).
    //
    // Window half-widths (same as Athena's egammaMiddleShape):
    //   E(7×7) = ±3.5 * deta_cell × ±3.5 * dphi_cell
    //   E(3×7) = ±1.5 * deta_cell × ±3.5 * dphi_cell
    //   E(3×3) = ±1.5 * deta_cell × ±1.5 * dphi_cell
    //   w_eta2 = ±1.5 * deta_cell × ±2.5 * dphi_cell
    // ======================================================================

    // Wrap delta-phi into [-pi, pi]
    inline double deltaPhi(double phi1, double phi2) {
        double d = phi1 - phi2;
        while (d >  M_PI) d -= 2 * M_PI;
        while (d < -M_PI) d += 2 * M_PI;
        return d;
    }

    // Sum cell energies within |eta-eta0| < halfEta and |phi-phi0| < halfPhi
    inline double EnergyPos(const std::vector<double>& E,
                            const std::vector<double>& eta,
                            const std::vector<double>& phi,
                            double eta0, double phi0,
                            double halfEta, double halfPhi) {
        double sum = 0;
        for (int k = 0; k < kClusterSize; ++k) {
            if (std::fabs(eta[k] - eta0) < halfEta &&
                std::fabs(deltaPhi(phi[k], phi0)) < halfPhi)
                sum += E[k];
        }
        return sum;
    }

    // Find the eta/phi of the hottest cell in the grid
    inline void findHottestCell(const std::vector<double>& E,
                                const std::vector<double>& eta,
                                const std::vector<double>& phi,
                                double& eta0, double& phi0) {
        double maxE = -1e30;
        eta0 = 0; phi0 = 0;
        for (int k = 0; k < kClusterSize; ++k) {
            if (E[k] > maxE) { maxE = E[k]; eta0 = eta[k]; phi0 = phi[k]; }
        }
    }

    inline double calcRetaPos(const std::vector<double>& E,
                              const std::vector<double>& eta,
                              const std::vector<double>& phi) {
        double eta0, phi0;
        findHottestCell(E, eta, phi, eta0, phi0);
        double E277 = EnergyPos(E, eta, phi, eta0, phi0,
                                3.5 * kDeltaEta, 3.5 * kDeltaPhi);
        double E237 = EnergyPos(E, eta, phi, eta0, phi0,
                                1.5 * kDeltaEta, 3.5 * kDeltaPhi);
        if (E277 <= 0) return -999.0;
        return E237 / E277;
    }

    inline double calcRphiPos(const std::vector<double>& E,
                              const std::vector<double>& eta,
                              const std::vector<double>& phi) {
        double eta0, phi0;
        findHottestCell(E, eta, phi, eta0, phi0);
        double E237 = EnergyPos(E, eta, phi, eta0, phi0,
                                1.5 * kDeltaEta, 3.5 * kDeltaPhi);
        double E233 = EnergyPos(E, eta, phi, eta0, phi0,
                                1.5 * kDeltaEta, 1.5 * kDeltaPhi);
        if (E237 <= 0) return -999.0;
        return E233 / E237;
    }

    inline double calcWeta2Pos(const std::vector<double>& E,
                               const std::vector<double>& cellEta,
                               const std::vector<double>& cellPhi,
                               float clusterEta, float etacell) {
        double eta0, phi0;
        findHottestCell(E, cellEta, cellPhi, eta0, phi0);
        double halfEta = 1.5 * kDeltaEta;
        double halfPhi = 2.5 * kDeltaPhi;
        double sumE = 0, sumEeta = 0, sumEeta2 = 0;
        for (int k = 0; k < kClusterSize; ++k) {
            if (std::fabs(cellEta[k] - eta0) < halfEta &&
                std::fabs(deltaPhi(cellPhi[k], phi0)) < halfPhi) {
                double e = E[k];
                double etaRel = cellEta[k] - eta0;
                sumE     += e;
                sumEeta  += e * etaRel;
                sumEeta2 += e * etaRel * etaRel;
            }
        }
        if (sumE <= 0) return -999.0;
        double var = sumEeta2 / sumE - std::pow(sumEeta / sumE, 2);
        if (var < 0) return -999.0;
        double raw = std::sqrt(var);
        return weta2Correct(clusterEta, etacell, raw);
    }

}  // namespace datamc

#endif // DATA_MC_CONFIG_H
