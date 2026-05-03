#ifndef DATA_MC_CONFIG_LAYER1_H
#define DATA_MC_CONFIG_LAYER1_H

///////////////////////////////////////////////////////////////////////////////
// config_layer1.h
//
// Layer-1 (EM strip) configuration for the cell-energy reweighting pipeline.
//
// Cell array from NTupleMaker:
//   photon.7x11ClusterLr1E     — 112 doubles (56 eta × 2 phi, row-major)
//   photon.7x11ClusterLr1Eta
//   photon.7x11ClusterLr1Phi
//   photon.7x11ClusterLr1Size  — number of non-padded cells
//
// Padded cells have (Eta=0, Phi=0, E=0).
//
// Analysis grid:
//   We use the FULL 56 × 2 = 112-cell cluster (same as xAOD shower shapes).
//   The hot (non-padded) strip is identified per event and used as the
//   reference index for strip-based widths and peak finding. This matches
//   the philosophy of the Layer-2 pipeline (correct cell energies on the
//   full cluster, recompute shower shapes from the same cluster) and is
//   required for w_s,tot (±20 strips) and ΔE/E_ratio (which scan the full
//   strip cluster for a second peak).
//
// Shower shapes computed from the 112-cell full grid:
//   weta1, w_s,tot, f_side, ΔE, E_ratio
//
// Selection, event weights, binning and scenario logic are inherited from
// config.h (Layer 2 configuration). Only Lr2-specific cluster-health is
// replaced by isHealthyL1Cluster; the Lr2 size cut inside passSelection is
// kept as a generic photon-quality cut.
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace config_l1 {

    using namespace config;  // Selection, binning, passSelection, MC weights

    // ======================================================================
    // Layer-1 cluster geometry
    // ======================================================================
    const int kLrEtaSize     = 56;
    const int kLrPhiSize     = 2;
    const int kLrClusterSize = kLrEtaSize * kLrPhiSize;   // 112

    // Full grid convention: k = iEta * kLrPhiSize + iPhi
    //   iEta ∈ [0, 55], iPhi ∈ {0, 1}
    // The hot strip index (hotEtaIdx) is found per event.
    const int kL1EtaSize  = kLrEtaSize;
    const int kL1PhiSize  = kLrPhiSize;
    const int kL1GridSize = kLrClusterSize;                // 112

    // Half-window for w_s,tot in L2 eta-coordinate units (5 L2 cells × 0.025 / 4).
    // Matches the deta*neta/4 = 0.025*5/4 = 0.03125 window used by ATLAS setWstot.
    const double kWsTotEtaHalf = 0.03125;

    // Relative grid used to accumulate cell-fraction profiles and per-cell
    // corrections, indexed by eta_rel = iEta - hotEtaIdx ∈ [-55, +55].
    // Size 111 × 2 = 222 covers every (hotEtaIdx, iEta) pair so no event
    // needs to be rejected for geometric reasons.
    //   kRel = (eta_rel + kRelEtaHalf) * kRelPhiSize + phi_rel
    //   hot strip is at eta_rel = 0 (kRelCenterEta = kRelEtaHalf).
    const int kRelEtaHalf    = kL1EtaSize - 1;             // 55
    const int kRelEtaSize    = 2 * kRelEtaHalf + 1;        // 111
    const int kRelPhiSize    = kL1PhiSize;                 // 2
    const int kRelGridSize   = kRelEtaSize * kRelPhiSize;  // 222
    const int kRelCenterEta  = kRelEtaHalf;                // 55

    // Peripheral-cell threshold for M2 regularisation
    const double kMinFracForM2_L1 = 0.002;

    // ======================================================================
    // Layer-1 branch names
    // ======================================================================
    const char* const kL1CellBranch     = "photon.7x11ClusterLr1E";
    const char* const kL1CellEtaBranch  = "photon.7x11ClusterLr1Eta";
    const char* const kL1CellPhiBranch  = "photon.7x11ClusterLr1Phi";
    const char* const kL1CellSizeBranch = "photon.7x11ClusterLr1Size";

    // Stored shower shapes (fudged for MC, standard for data)
    const char* const kWeta1Branch   = "photon.weta1";
    const char* const kWstot1Branch  = "photon.wtots1";
    const char* const kFsideBranch   = "photon.fracs1";
    const char* const kDeltaEBranch  = "photon.deltae";
    const char* const kEratioBranch  = "photon.eratio";

    // Unfudged (MC only)
    const char* const kWeta1UnfBranch   = "photon.unfudged_weta1";
    const char* const kWstot1UnfBranch  = "photon.unfudged_wtots1";
    const char* const kFsideUnfBranch   = "photon.unfudged_fracs1";
    const char* const kDeltaEUnfBranch  = "photon.unfudged_deltae";
    const char* const kEratioUnfBranch  = "photon.unfudged_eratio";

    // ======================================================================
    // Cell utilities
    // ======================================================================

    // A padded cell has (eta=0, phi=0, E=0).
    inline bool isPadded(double eta, double phi) {
        return eta == 0.0 && phi == 0.0;
    }

    // Build the full 56×2 grid (all 112 cells, padded cells kept as 0 energy)
    // and return the index of the hottest non-padded strip (hotEtaIdx in
    // [0,55]). The return value is false only for gross structural problems
    // (wrong array size, all cells padded, zero total energy).
    //
    // full[112] stores energies row-major: k = iEta * kL1PhiSize + iPhi.
    inline bool buildFullGrid(const std::vector<double>& cellE,
                              const std::vector<double>& cellEta,
                              const std::vector<double>& cellPhi,
                              std::vector<double>& full,
                              int& hotEtaIdx, int& hotPhiIdx) {
        if ((int)cellE.size()   != kL1GridSize) return false;
        if ((int)cellEta.size() != kL1GridSize) return false;
        if ((int)cellPhi.size() != kL1GridSize) return false;

        full.assign(kL1GridSize, 0.0);

        int kHot = -1;
        double eHot = -1e30;
        for (int k = 0; k < kL1GridSize; ++k) {
            // Keep all cell energies (including padded ones, which are 0).
            full[k] = cellE[k];
            if (isPadded(cellEta[k], cellPhi[k])) continue;
            if (cellE[k] > eHot) { eHot = cellE[k]; kHot = k; }
        }
        if (kHot < 0 || eHot <= 0) return false;

        hotEtaIdx = kHot / kL1PhiSize;
        hotPhiIdx = kHot % kL1PhiSize;
        return true;
    }

    inline bool isHealthyL1Cluster(const std::vector<double>& cellE,
                                   const std::vector<double>& cellEta,
                                   const std::vector<double>& cellPhi) {
        std::vector<double> full;
        int he = -1, hp = -1;
        return buildFullGrid(cellE, cellEta, cellPhi, full, he, hp);
    }

    // Returns true if the (iEta, iPhi) cell is a real (non-padded) cell.
    inline bool isRealCell(const std::vector<double>& cellEta,
                           const std::vector<double>& cellPhi,
                           int iEta, int iPhi) {
        int k = iEta * kL1PhiSize + iPhi;
        return !isPadded(cellEta[k], cellPhi[k]);
    }

    // ======================================================================
    // ATLAS-matching geometric helpers for shower shape corrections
    // ======================================================================

    inline double stripGranularity(double aeta) {
        if (aeta < 1.8) return 0.025 / 8.0;
        if (aeta < 2.0) return 0.025 / 6.0;
        if (aeta < 2.4) return 0.025 / 4.0;
        if (aeta < 2.5) return 0.025 / 1.0;
        return 0.025 / 8.0;
    }

    inline double relPosition1(double eta, double etacell) {
        double aeta   = std::fabs(eta);
        double dgra   = stripGranularity(aeta);
        double etapos = std::fabs(eta - etacell - dgra / 2.0);
        return std::fmod(etapos, dgra) / dgra - 0.5;
    }

    inline double correctWs3(double eta, double etacell, double ws3) {
        if (ws3 < 0) return ws3;
        double aeta = std::fabs(eta);
        double u  = relPosition1(eta, etacell);
        double u2 = u * u;
        double u4 = u2 * u2;
        if (aeta < 1.0)  return ws3 - 0.76 * u2;
        if (aeta < 1.45) return ws3 - 0.85 * u2 + 1.9 * u4;
        if (aeta < 1.5)  return ws3;
        if (aeta < 1.8)  return ws3 - 0.85 * u2 + 1.9 * u4;
        if (aeta < 2.0)  return ws3 - 0.84 * u2;
        if (aeta < 2.5)  return ws3 - 0.40 * u2 - 2.1 * u4;
        return ws3;
    }

    // ======================================================================
    // Shower-shape calculators on the full 56×2 grid
    //
    // Index convention for full[]:  k = iEta * kL1PhiSize + iPhi
    //   iEta ∈ [0, 55], iPhi ∈ {0, 1}
    // The hot strip is at iEta = hotEtaIdx.
    //
    // Widths use dimensionless strip indices (i - hotEtaIdx), matching xAOD.
    // ΔE and E_ratio use absolute cell energies.
    // ======================================================================

    // Sum over iEta ∈ [iLo, iHi] (clamped to [0,55]), both phi.
    inline double sumEAbs(const std::vector<double>& full, int iLo, int iHi) {
        if (iLo < 0) iLo = 0;
        if (iHi > kL1EtaSize - 1) iHi = kL1EtaSize - 1;
        double s = 0;
        for (int i = iLo; i <= iHi; ++i) {
            int base = i * kL1PhiSize;
            s += full[base] + full[base + 1];
        }
        return s;
    }

    // 1D eta projection (sum over 2 phi cells) of the full cluster.
    inline std::vector<double> projectEtaFull(const std::vector<double>& full) {
        std::vector<double> p(kL1EtaSize, 0.0);
        for (int i = 0; i < kL1EtaSize; ++i)
            p[i] = full[i * kL1PhiSize] + full[i * kL1PhiSize + 1];
        return p;
    }

    // w_eta_1: energy-weighted strip-index width over 3 strips × 2 phi
    // centered on the hot strip, with ATLAS position-within-strip correction.
    // Matches ATLAS setWs3: width = sqrt(sum(E*d^2) / sum(E))  — NO mean^2
    // subtraction (i.e. width is taken about the hot-strip index, not about
    // the energy-weighted mean).
    inline double calcWeta1(const std::vector<double>& full, int hotEtaIdx,
                            const std::vector<double>& cellEta,
                            double clusterEta1, double hotEta1) {
        double sw = 0, s2 = 0;
        for (int d = -1; d <= 1; ++d) {
            int i = hotEtaIdx + d;
            if (i < 0 || i >= kL1EtaSize) continue;
            int base = i * kL1PhiSize;
            double e = full[base] + full[base + 1];
            if (e <= 0) continue;
            sw += e; s2 += e * d * d;
        }
        if (sw <= 0) return -999.0;
        double var = s2 / sw;
        if (var < 0) return -999.0;
        double ws3 = std::sqrt(var);
        return correctWs3(clusterEta1, hotEta1, ws3);
    }

    // w_s,tot: energy-weighted strip-index width over all strips within
    // ±kWsTotEtaHalf (0.03125) in eta from the L2 hottest cell.
    // Matches ATLAS setWstot (deta*neta/4 = 0.025*5/4 window).
    inline double calcWsTot(const std::vector<double>& full, int hotEtaIdx,
                            const std::vector<double>& cellEta, double etaMaxL2) {
        double sw = 0, s2 = 0;
        for (int i = 0; i < kL1EtaSize; ++i) {
            double eta_i = cellEta[i * kL1PhiSize];
            if (eta_i == 0.0) continue;
            if (std::fabs(eta_i - etaMaxL2) >= kWsTotEtaHalf) continue;
            int base = i * kL1PhiSize;
            double e = full[base] + full[base + 1];
            if (e <= 0) continue;
            double d = static_cast<double>(i - hotEtaIdx);
            sw += e;
            s2 += e * d * d;
        }
        if (sw <= 0) return -999.0;
        double var = s2 / sw;
        if (var < 0) return -999.0;
        return std::sqrt(var);
    }

    // f_side: ATLAS setFside.
    //   fside = (sum_{7-strip} E/gracell) / (E_hot/g_hot + E_left/g_left
    //           + E_right/g_right) - 1
    // with the M.S. 60 MeV threshold on the hot-strip energy (below threshold
    // → fside = 0, matching ATLAS).
    inline double calcFside(const std::vector<double>& full, int hotEtaIdx,
                            const std::vector<double>& cellEta) {
        std::vector<double> p = projectEtaFull(full);
        int N = (int)p.size();
        auto gc = [&](int i) -> double {
            return stripGranularity(std::fabs(cellEta[i * kL1PhiSize]));
        };
        double e1 = p[hotEtaIdx];
        if (e1 <= 0) return -999.0;
        // ATLAS M.S. threshold (60 MeV) on hot-strip energy; below it,
        // fside is left at its default 0.
        const double kFsideThreshold = 60.0;  // MeV (xAOD energies are MeV)
        if (e1 <= kFsideThreshold) return 0.0;
        int ileft = hotEtaIdx - 1;
        while (ileft >= 0  && p[ileft] == 0) --ileft;
        int iright = hotEtaIdx + 1;
        while (iright < N && p[iright] == 0) ++iright;
        double eleft  = (ileft  >= 0) ? p[ileft]  : 0.0;
        double eright = (iright <  N) ? p[iright] : 0.0;
        double fracm = 0.0;
        int nlo = std::max(0, hotEtaIdx - 3);
        int nhi = std::min(N - 1, hotEtaIdx + 3);
        for (int i = nlo; i <= nhi; ++i) {
            double g = gc(i);
            if (g > 0) fracm += p[i] / g;
        }
        double denom = e1 / gc(hotEtaIdx);
        if (ileft  >= 0) denom += eleft  / gc(ileft);
        if (iright <  N) denom += eright / gc(iright);
        if (denom <= 0) return -999.0;
        return fracm / denom - 1.0;
    }

    // Peak finder on the 1D eta profile of the full cluster.
    // Iterates over all 56 strips. Strips with zero energy (typically padded
    // cells at the edges) are naturally excluded since they cannot be maxima.
    struct PeakInfo {
        int    iMax1 = -1;
        int    iMax2 = -1;
        int    iMin  = -1;
        double eMax1 = 0.0;
        double eMax2 = 0.0;
        double eMin  = 0.0;
        bool   ok    = false;
        bool   hasSecondary = false;
    };

    inline PeakInfo findPeaks(const std::vector<double>& full,
                              const std::vector<double>& cellEta) {
        PeakInfo info;
        std::vector<double> p = projectEtaFull(full);
        int N = (int)p.size();  // 56

        // Energy-density array: d[i] = p[i] / gracell[i]. Used for ranking
        // local maxima and finding the inter-peak minimum (matches ATLAS
        // setEmax2 / setEmin which compare e/gracell).
        std::vector<double> d(N, 0.0);
        for (int i = 0; i < N; ++i) {
            double g = stripGranularity(std::fabs(cellEta[i * kL1PhiSize]));
            if (g > 0) d[i] = p[i] / g;
        }

        // Primary maximum: hottest strip in raw energy (ATLAS setEmax).
        info.iMax1 = 0; info.eMax1 = p[0];
        for (int i = 1; i < N; ++i)
            if (p[i] > info.eMax1) { info.eMax1 = p[i]; info.iMax1 = i; }
        if (info.eMax1 <= 0) return info;

        // Secondary maximum (esec1): largest local max (by density) at
        // least 2 strips away, stored as raw 1-strip energy.
        double dMax2 = -1e30;
        for (int i = 1; i < N - 1; ++i) {
            if (std::abs(i - info.iMax1) < 2) continue;
            if (d[i] <= 0) continue;
            if (d[i] <= d[i - 1] || d[i] <= d[i + 1]) continue;
            if (d[i] > dMax2) { dMax2 = d[i]; info.iMax2 = i; }
        }
        if (info.iMax2 < 0) { info.ok = true; info.eMax2 = 0; return info; }
        info.eMax2 = p[info.iMax2];

        // Minimum strictly between the two peaks (ATLAS setEmin uses
        // [min+1, max-1]). If the peaks are adjacent the range is empty
        // and emins1 is 0.
        int lo = std::min(info.iMax1, info.iMax2) + 1;
        int hi = std::max(info.iMax1, info.iMax2) - 1;
        info.iMin = -1;
        info.eMin = 0.0;
        if (lo <= hi) {
            double dMin = std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i) {
                if (d[i] < dMin) { dMin = d[i]; info.iMin = i; }
            }
            if (info.iMin >= 0) info.eMin = p[info.iMin];
        }
        info.ok = true; info.hasSecondary = true;
        return info;
    }

    // ΔE = esec1 - emins1 (MeV) — the 1-strip 2nd-max energy minus the
    // 1-strip minimum energy between the two peaks (ATLAS
    // EMShowerBuilder::FillEMShowerShape stores DeltaE = esec1 - emins1).
    // Convention: if no secondary maximum, 0.
    inline double calcDeltaE(const std::vector<double>& full,
                             const std::vector<double>& cellEta) {
        PeakInfo info = findPeaks(full, cellEta);
        if (!info.ok) return -999.0;
        if (!info.hasSecondary) return 0.0;
        return info.eMax2 - info.eMin;
    }

    // E_ratio = (E_max1 - E_max2) / (E_max1 + E_max2).
    // Convention: if no secondary maximum, 1.
    inline double calcEratio(const std::vector<double>& full,
                             const std::vector<double>& cellEta) {
        PeakInfo info = findPeaks(full, cellEta);
        if (!info.ok) return -999.0;
        if (!info.hasSecondary) return 1.0;
        double denom = info.eMax1 + info.eMax2;
        if (denom <= 0) return -999.0;
        return (info.eMax1 - info.eMax2) / denom;
    }

    // ======================================================================
    // Histogram ranges for each variable (for booking TH1D)
    // ======================================================================
    struct SSRange { double lo, hi; };
    inline SSRange rangeWeta1 () { return {0.30, 0.90}; }  // strip units
    inline SSRange rangeWstot () { return {0.00, 8.00}; }  // strip units
    inline SSRange rangeFside () { return {0.00, 1.00}; }
    inline SSRange rangeDeltaE() { return {0.0, 500.0}; }  // MeV
    inline SSRange rangeEratio() { return {0.50, 1.05}; }

    // ======================================================================
    // M3 / M4 helpers
    //
    // M3 = per-cell quantile transport. For each (eta, pT, cell) we hold an
    //      empirical CDF from data and from MC, and map each MC f_k through
    //      f_k' = F_data^{-1}(F_mc(f_k)).
    //
    // M4 = M3 followed by a rank-2 position reshuffle: the current rank-2
    //      cell's energy is swapped with a cell drawn from the data
    //      distribution of rank-2 positions (relative to the hot strip).
    //
    // Both methods inherit the (eta, pT) binning — more bins = more
    // kinematic fidelity but fewer events per (eta, pT, cell) triple, so
    // CDF tails degrade.  Fallback to identity when statistics are too low.
    // ======================================================================

    // Minimum events per TH1D before trusting its empirical CDF.
    const double kMinCDFCounts = 100.0;

    // Quantile transport for a single value f using the MC source CDF and
    // the data target CDF.  Returns the input unchanged when either
    // histogram has too few entries (numerically unstable tails).
    //
    // Implementation notes:
    //   - Uses TH1::GetIntegral() which returns the (normalised to 1) CDF
    //     evaluated at the right edge of each bin, plus an overflow term.
    //   - Bilinear interpolation is avoided: we clip to (qLo, qHi) = (0.5%,
    //     99.5%) to tame tail instability.
    //   - If f falls outside the MC histogram range, return f unchanged.
    inline double applyQuantileTransport(TH1D* hMC, TH1D* hData, double f) {
        if (!hMC || !hData) return f;
        if (hMC->GetEntries()  < kMinCDFCounts) return f;
        if (hData->GetEntries() < kMinCDFCounts) return f;

        // Source CDF at f.
        int bin = hMC->FindFixBin(f);
        if (bin < 1 || bin > hMC->GetNbinsX()) return f;
        double sumTot = hMC->Integral(1, hMC->GetNbinsX());
        if (sumTot <= 0) return f;
        double sumBelow = hMC->Integral(1, bin - 1);
        double bc = hMC->GetBinContent(bin);
        double lo = hMC->GetBinLowEdge(bin);
        double w  = hMC->GetBinWidth(bin);
        double frac = (w > 0) ? (f - lo) / w : 0.0;
        double u = (sumBelow + frac * bc) / sumTot;

        // Clip to avoid extrapolating the inverse CDF into empty tails.
        const double qLo = 0.005, qHi = 0.995;
        if (u < qLo) u = qLo;
        if (u > qHi) u = qHi;

        // Target inverse-CDF via GetQuantiles (linear interpolation between
        // bin edges).  GetQuantiles needs a non-const pointer, so const_cast.
        double pArr[1] = {u};
        double qArr[1] = {0.0};
        hData->GetQuantiles(1, qArr, pArr);
        return qArr[0];
    }

    // Find the absolute indices of the top-2 energy cells on the full 56x2
    // grid.  Returns false if fewer than 2 non-padded cells carry positive
    // energy.  Conventions:
    //   k = iEta * kL1PhiSize + iPhi
    //   kMax1 = hot cell (guaranteed by buildFullGrid)
    //   kMax2 = largest energy cell excluding kMax1
    inline bool findTop2Abs(const std::vector<double>& full,
                            int hotEtaIdx, int hotPhiIdx,
                            int& kMax1, int& kMax2) {
        kMax1 = hotEtaIdx * kL1PhiSize + hotPhiIdx;
        kMax2 = -1;
        double e2 = -1e30;
        for (int k = 0; k < kL1GridSize; ++k) {
            if (k == kMax1) continue;
            if (full[k] > e2) { e2 = full[k]; kMax2 = k; }
        }
        return (kMax2 >= 0 && e2 > 0);
    }

}  // namespace config_l1

#endif  // DATA_MC_CONFIG_LAYER1_H
