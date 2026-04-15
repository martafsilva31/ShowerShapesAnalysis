#ifndef DATA_MC_CONFIG_H
#define DATA_MC_CONFIG_H

///////////////////////////////////////////////////////////////////////////////
// config.h
//
// Shared configuration for the cell-energy reweighting pipeline.
//
// Two correction methods:
//   M1  — Flat shift:          E'_k = E_k + delta_k * E_total
//   M2  — Shift + stretch:     E'_k = E_total * shift_k + stretch_k * E_k
//          (equivalent to Francisco's "method1" in photoncellbasedrw)
//
// Data files:   egam3.root (Z->eeg),  egam4.root (Z->mumug)
// MC files:     mc_eegamma.root,      mc_mumugamma.root
//               use channel "llgamma" to chain both simultaneously
// Location:     /dcache/atlas/mfernand/qt_ntuples/data24/
//
// Reference: ATL-COM-PHYS-2025-662, ATL-COM-PHYS-2021-640
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <map>

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

namespace config {

    // ======================================================================
    // Cluster geometry  (EM Layer 2, 7 eta x 11 phi = 77 cells)
    // ======================================================================
    const int kEtaSize     = 7;
    const int kPhiSize     = 11;
    const int kClusterSize = kEtaSize * kPhiSize;   // 77
    const int kCentralCell = 3 * kPhiSize + 5;       // 38
    const double kDeltaEta = 0.025;
    const double kDeltaPhi = 0.0245;

    // Peripheral-cell threshold for M2 regularisation:
    // cells with <f_k> below this fall back to M1 (stretch=1).
    const double kMinFracForM2 = 0.002;

    // ======================================================================
    // Eta binning (14 bins, crack at bin 8)
    // ======================================================================
    const int kNEtaBins = 14;
    const double kEtaLimits[kNEtaBins + 1] = {
        0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.37,
        1.52, 1.6, 1.80, 2.0, 2.2, 2.4
    };

    // ======================================================================
    // pT binning (6 bins, for future use)
    // ======================================================================
    const int kNPtBins = 6;
    const double kPtLimits[kNPtBins + 1] = {
        10, 15, 20, 25, 30, 40, 1000   // GeV
    };

    // ======================================================================
    // Selection scenario
    //
    // Baseline = Francisco's defaults:
    //   pT > 10 GeV,  noID,  loose isolation,  no lepton-lepton DR cut
    // ======================================================================
    struct Selection {
        double photonPtMin     = 10.0;      // GeV
        double etaMax          = 2.37;
        double crackEtaMin     = 1.37;
        double crackEtaMax     = 1.52;
        double mllMin          = 40.0;      // GeV
        double mllMax          = 83.0;      // GeV
        double mllgMin         = 80.0;      // GeV
        double mllgMax         = 100.0;     // GeV
        double dRPhotonLepMin  = 0.4;
        double dRLeptonLepMin  = 0.0;       // Francisco doesn't cut on this
        bool   requireTightID  = false;     // Francisco: noID
        bool   requireLooseIso = true;      // Francisco: loose
        bool   requireTightIso = false;
        bool   unconvertedOnly = true;
        bool   convertedOnly   = false;
        bool   requireTruthMC  = false;     // Francisco: no truth cut
    };

    inline Selection baselineSelection() {
        return Selection{};
    }

    inline Selection tightIDSelection() {
        Selection s;
        s.requireTightID = true;
        return s;
    }

    inline Selection isoTightSelection() {
        Selection s;
        s.requireTightIso = true;
        s.requireLooseIso = false;
        return s;
    }

    inline Selection noIsoSelection() {
        Selection s;
        s.requireLooseIso = false;
        return s;
    }

    inline Selection tightIDTightIsoSelection() {
        Selection s;
        s.requireTightID  = true;
        s.requireTightIso = true;
        s.requireLooseIso = false;
        return s;
    }

    inline Selection convertedSelection() {
        Selection s;
        s.unconvertedOnly = false;
        s.convertedOnly   = true;
        return s;
    }

    inline Selection allConvSelection() {
        Selection s;
        s.unconvertedOnly = false;
        s.convertedOnly   = false;
        return s;
    }

    inline Selection isoTightConvertedSelection() {
        Selection s;
        s.requireTightIso = true;
        s.requireLooseIso = false;
        s.unconvertedOnly = false;
        s.convertedOnly   = true;
        return s;
    }

    inline Selection isoTightAllConvSelection() {
        Selection s;
        s.requireTightIso = true;
        s.requireLooseIso = false;
        s.unconvertedOnly = false;
        s.convertedOnly   = false;
        return s;
    }

    inline Selection getSelection(const std::string& name) {
        if (name == "baseline")            return baselineSelection();
        if (name == "converted")           return convertedSelection();
        if (name == "all_conv")            return allConvSelection();
        if (name == "tight_id")            return tightIDSelection();
        if (name == "iso_tight")           return isoTightSelection();
        if (name == "no_iso")              return noIsoSelection();
        if (name == "tight_id_tight_iso")  return tightIDTightIsoSelection();
        if (name == "iso_tight_converted")   return isoTightConvertedSelection();
        if (name == "iso_tight_all_conv")    return isoTightAllConvSelection();
        std::cerr << "ERROR: Unknown scenario '" << name << "'" << std::endl;
        return baselineSelection();
    }

    // ======================================================================
    // Human-readable labels for plots
    // ======================================================================
    inline const char* channelLabel(const char* channel) {
        TString ch(channel);
        if (ch == "eegamma")    return "Z#rightarrowee#gamma";
        if (ch == "mumugamma")  return "Z#rightarrow#mu#mu#gamma";
        if (ch == "llgamma")    return "Z#rightarrowll#gamma";
        return channel;
    }

    inline const char* scenarioLabel(const char* scenario) {
        TString sc(scenario);
        if (sc == "baseline")    return "Unconverted #gamma, no ID, loose iso";
        if (sc == "converted")   return "Converted #gamma, no ID, loose iso";
        if (sc == "all_conv")    return "All #gamma, no ID, loose iso";
        if (sc == "tight_id")    return "Unconverted #gamma, tight ID, loose iso";
        if (sc == "iso_tight")             return "Unconverted #gamma, no ID, tight iso";
        if (sc == "iso_tight_converted")   return "Converted #gamma, no ID, tight iso";
        if (sc == "iso_tight_all_conv")    return "All #gamma, no ID, tight iso";
        if (sc == "no_iso")                return "Unconverted #gamma, no ID, no iso";
        return scenario;
    }

    // ======================================================================
    // Branch names for data24 / MC23e ntuples (from NTupleMaker)
    // ======================================================================
    const char* const kTreeName = "tree";

    // Cell energies and positions (Layer 2, 7x11 cluster)
    const char* const kCellBranch     = "photon.7x11ClusterLr2E";
    const char* const kCellEtaBranch  = "photon.7x11ClusterLr2Eta";
    const char* const kCellPhiBranch  = "photon.7x11ClusterLr2Phi";
    const char* const kCellSizeBranch = "photon.7x11ClusterLr2Size";

    // Cluster-level kinematics
    const char* const kClusterEta2Branch    = "photon_cluster.eta2";
    const char* const kClusterEtamax2Branch = "photon_cluster.etamax2";
    const char* const kPhotonPtBranch       = "photon.pt";
    const char* const kPhotonEtaBranch      = "photon.eta";
    const char* const kPhotonPhiBranch      = "photon.phi";
    const char* const kPhotonEBranch        = "photon_cluster.e";

    // Mass variables
    const char* const kMllBranch  = "mll";
    const char* const kMllgBranch = "mllg";

    // Photon properties
    const char* const kIsConvBranch  = "photon.isconv";
    const char* const kTightIDBranch = "photon.tightOffInc";
    const char* const kLooseIDBranch = "photon.looseOff";
    const char* const kIsoLooseBranch = "photon.isoloose";
    const char* const kIsoTightBranch = "photon.isotight";

    // Lepton-photon / lepton-lepton separation
    const char* const kDRLepton1Branch = "lepton1.dr_ph";
    const char* const kDRLepton2Branch = "lepton2.dr_ph";
    const char* const kDRllBranch      = "dRll";

    // MC weights
    const char* const kMCWeightBranch  = "event.mcwgt";
    const char* const kPUWeightBranch  = "event.muwgt";
    const char* const kXSecBranch      = "event.xsec";
    const char* const kBSWeightBranch  = "event.beamspotwgt";
    const char* const kLepton1SFBranch = "lepton1.SF";
    const char* const kLepton2SFBranch = "lepton2.SF";
    const char* const kPhotonIDSFBranch  = "photon.ID_sf";
    const char* const kPhotonIsoSFBranch = "photon.sf";   // isolation SF from NTupleMaker acc_sf("sf")

    // MC truth
    const char* const kTruthMatchBranch = "photon.istruthmatch";

    // Stored shower shapes (fudged for MC, standard for data)
    const char* const kRetaBranch  = "photon.reta";
    const char* const kRphiBranch  = "photon.rphi";
    const char* const kWeta2Branch = "photon.weta2";

    // Unfudged shower shapes (MC only)
    const char* const kRetaUnfBranch  = "photon.unfudged_reta";
    const char* const kRphiUnfBranch  = "photon.unfudged_rphi";
    const char* const kWeta2UnfBranch = "photon.unfudged_weta2";

    // ======================================================================
    // Inline helpers
    // ======================================================================

    inline int findEtaBin(double abseta) {
        for (int n = 0; n < kNEtaBins; ++n)
            if (abseta >= kEtaLimits[n] && abseta < kEtaLimits[n + 1])
                return n;
        return -1;
    }

    inline int findPtBin(double pt) {
        for (int n = 0; n < kNPtBins; ++n)
            if (pt >= kPtLimits[n] && pt < kPtLimits[n + 1])
                return n;
        return -1;
    }

    // ======================================================================
    // Combined event selection
    // ======================================================================
    inline bool passSelection(const Selection& sel,
                              double photonPt,   // GeV
                              double clusterEta,
                              double mll,        // GeV
                              double mllg,       // GeV
                              double dR_l1_ph,
                              double dR_l2_ph,
                              double dR_ll,
                              int    cellSize,
                              bool   isConv,
                              int    isTight,
                              int    isLooseIso,
                              int    isTightIso,
                              bool   isTruthMatch,
                              bool   isMC) {
        if (sel.unconvertedOnly && isConv) return false;
        if (sel.convertedOnly && !isConv) return false;
        if (photonPt < sel.photonPtMin) return false;
        double abseta = std::fabs(clusterEta);
        if (abseta > sel.etaMax) return false;
        if (abseta > sel.crackEtaMin && abseta < sel.crackEtaMax) return false;
        if (mll < sel.mllMin || mll > sel.mllMax) return false;
        if (mllg < sel.mllgMin || mllg > sel.mllgMax) return false;
        if (dR_l1_ph < sel.dRPhotonLepMin) return false;
        if (dR_l2_ph < sel.dRPhotonLepMin) return false;
        if (sel.dRLeptonLepMin > 0 && dR_ll < sel.dRLeptonLepMin) return false;
        if (cellSize != kClusterSize) return false;
        if (sel.requireTightID && isTight != 1) return false;
        if (sel.requireLooseIso && isLooseIso != 1) return false;
        if (sel.requireTightIso && isTightIso != 1) return false;
        if (sel.requireTruthMC && isMC && !isTruthMatch) return false;
        return true;
    }

    // ======================================================================
    // Cluster quality
    // ======================================================================
    inline bool isHealthyCluster(const std::vector<double>& cells) {
        if ((int)cells.size() != kClusterSize) return false;
        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k)
            Etot += cells[k];
        if (Etot <= 0) return false;
        double centralE = cells[kCentralCell];
        for (int k = 0; k < kClusterSize; ++k) {
            if (k == kCentralCell) continue;
            if (cells[k] > centralE) return false;
        }
        return true;
    }

    // ======================================================================
    // TChain helper (comma-separated file list)
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

    // ======================================================================
    // Shower shape calculations — index-based windows
    // ======================================================================
    inline double Energy(int eta, int phi,
                         const std::vector<double>& clusterEnergy) {
        int etaMin = kEtaSize - (kEtaSize + eta) / 2;
        int etaMax = kEtaSize - (kEtaSize - eta) / 2;
        int phiMin = kPhiSize - (kPhiSize + phi) / 2;
        int phiMax = kPhiSize - (kPhiSize - phi) / 2;
        double sumE = 0;
        for (int e = etaMin; e < etaMax; ++e)
            for (int p = phiMin; p < phiMax; ++p)
                sumE += clusterEnergy[p + kPhiSize * e];
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
    // ======================================================================
    inline double weta2RelPosition(float eta, float etacell) {
        double x = std::fabs(eta - etacell - 0.025 / 2.);
        return std::fmod(x, 0.025) / 0.025;
    }

    inline float weta2Correct(float eta, float etacell, float weta2_raw) {
        if (eta == -999.f || etacell == -999.f) return weta2_raw;
        float aeta = std::fabs(eta);
        float u = weta2RelPosition(eta, etacell);

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

    inline double calcWeta2(const std::vector<double>& cells,
                            float eta, float etacell) {
        double raw = calcWeta2Raw(cells);
        if (raw < 0) return raw;
        return weta2Correct(eta, etacell, raw);
    }

    // ======================================================================
    // MC event weight — absolute normalisation
    //
    // Per-event weight = mcwgt * muwgt * (σ * L / sumW) * SFs
    //
    // Cross-sections from PMG: PMGxsecDB_mc23.txt
    //   700770  Sh_2214_eegamma     σ = 102.60 pb   genFilterEff = 1.0
    //   700771  Sh_2214_mumugamma   σ = 102.59 pb   genFilterEff = 1.0
    //
    // Luminosity: data24 = 108 fb⁻¹ = 108000 pb⁻¹
    // ======================================================================
    const double kLumi_pb = 108.0e3;   // 108 fb⁻¹ in pb⁻¹

    // Cross-section per DSID [pb]
    inline double xsecForDSID(int dsid) {
        switch (dsid) {
            case 700770: return 102.60;   // Sh_2214_eegamma
            case 700771: return 102.59;   // Sh_2214_mumugamma
            default:
                std::cerr << "WARNING: unknown DSID " << dsid
                          << ", using xsec=1 pb\n";
                return 1.0;
        }
    }

    // Sum-of-weights from AMI (EVNT level, full MC sample).
    //
    //   DSID 700770  Sh_2214_eegamma   sumW = 5.269047e+14  (totalEvents = 185,402,000)
    //   DSID 700771  Sh_2214_mumugamma sumW = 5.421893e+14  (totalEvents = 189,985,000)
    //
    // The mcFiles argument is kept for API compatibility but is not used.
    inline std::map<int, double> computeSumWeightsFromFiles(const char* /*mcFiles*/) {
        std::map<int, double> result;
        result[700770] = 5.269047e+14;   // Sh_2214_eegamma   (AMI EVNT)
        result[700771] = 5.421893e+14;   // Sh_2214_mumugamma (AMI EVNT)
        return result;
    }

    // Per-event normalisation factor: σ * L / sumW_DSID
    // Must be multiplied by mcwgt * muwgt * SFs
    inline double mcNormFactor(int dsid,
                               const std::map<int, double>& sumWmap) {
        auto it = sumWmap.find(dsid);
        if (it == sumWmap.end()) {
            std::cerr << "WARNING: no sumW for DSID " << dsid << "\n";
            return 0.0;
        }
        return xsecForDSID(dsid) * kLumi_pb / it->second;
    }

    inline double mcWeight(double mcwgt, double muwgt, double normFactor,
                           double l1sf, double l2sf, double phIsosf) {
        return mcwgt * muwgt * normFactor * l1sf * l2sf * phIsosf;
    }

}  // namespace config

#endif // DATA_MC_CONFIG_H
