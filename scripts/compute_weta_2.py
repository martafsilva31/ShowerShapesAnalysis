#!/usr/bin/env python3
"""
Compute w_eta_2 from Layer 2 cell energies and compare with stored values.

Computes weta2 using grid-index positions (3x5 window in centre of 7x11 grid),
then applies the egammaqweta2c position-dependent correction using
photon_cluster.etaSample2 and photon_cluster.etamax2 (matching Athena exactly).

Usage:
    python compute_weta_2.py -i ntuples/mc23e/mc23e_700770_Zeeg.root -o output/weta2_Zeeg.root
"""

import ROOT
import numpy as np
from utils import make_parser, extract_sample_label

# ---------------------------------------------------------------------------
# Compiled C++ computation functions
# ---------------------------------------------------------------------------
ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>

// --- Raw weta2 using actual stored cell eta positions (3x5 window in 7x11 grid) ---
// Matches Athena's egammaMiddleShape: uses cell->eta() (geometric center) for each cell.
double compute_weta2_raw(const std::vector<double>& cells, const std::vector<double>& etas) {
    if (cells.size() != 77 || etas.size() != 77) return -999.0;
    const int nphi = 11;
    double sumE = 0, sumEeta = 0, sumEeta2 = 0;
    for (int ie = 2; ie < 5; ie++) {
        for (int ip = 3; ip < 8; ip++) {
            int idx = ie * nphi + ip;
            double e   = cells[idx];
            double eta = etas[idx];
            sumE     += e;
            sumEeta  += e * eta;
            sumEeta2 += e * eta * eta;
        }
    }
    if (sumE <= 0) return -999.0;
    double var = sumEeta2 / sumE - std::pow(sumEeta / sumE, 2);
    return (var >= 0) ? std::sqrt(var) : -999.0;
}

// --- egammaqweta2c correction (exact copy from Athena 25.0.40) ---
namespace qweta2c {
    constexpr double P0A[3] = {0.0045,   0.005375, -0.0562};
    constexpr double P1A[3] = {-0.0016, -0.0215,    0.114 };
    constexpr double P2A[3] = {-0.0866,  0.0215,   -0.053 };

    constexpr double P0B[3] = {0.0039,   0.005075, -0.0324};
    constexpr double P1B[3] = {0.00816, -0.0203,    0.0653};
    constexpr double P2B[3] = {-0.145,   0.0203,   -0.0286};

    constexpr double P0C[3] = {0.0047,  0.0035,  0.0};
    constexpr double P1C[3] = {-0.0184, -0.0139, 0.0};
    constexpr double P2C[3] = {0.0180,  0.0137,  0.0};

    double RelPosition(double eta, double etacell) {
        if (eta == -999.) return -999;
        const double x = std::fabs(eta - etacell - 0.025 / 2.);
        return std::fmod(x, 0.025) / 0.025;
    }

    double Correct(double eta, double etacell, double weta2) {
        double aeta = std::fabs(eta);
        double etarel = RelPosition(eta, etacell);
        const double *P0, *P1, *P2;
        int idx;

        if (aeta < 0.8) {
            P0 = P0A; P1 = P1A; P2 = P2A;
        } else if (aeta < 1.8) {
            P0 = P0B; P1 = P1B; P2 = P2B;
        } else if (aeta < 2.0) {
            return weta2 - (P0C[0] + P1C[0]*etarel + P2C[0]*etarel*etarel);
        } else if (aeta < 2.5) {
            return weta2 - (P0C[1] + P1C[1]*etarel + P2C[1]*etarel*etarel);
        } else {
            return weta2;
        }

        if (etarel < 0.1) idx = 0;
        else if (etarel < 0.9) idx = 1;
        else idx = 2;
        return weta2 - (P0[idx] + P1[idx]*etarel + P2[idx]*etarel*etarel);
    }
}

// --- Compute weta2 with Athena-exact correction ---
double compute_weta2_corrected(const std::vector<double>& cells,
                               const std::vector<double>& etas,
                               double etaSample2, double etamax2) {
    double raw = compute_weta2_raw(cells, etas);
    if (raw < 0) return -999.0;
    return qweta2c::Correct(etaSample2, etamax2, raw);
}
""")


def main():
    parser = make_parser('Compute w_eta_2 from cell energies')
    args = parser.parse_args()

    sample = extract_sample_label(args.input)

    f = ROOT.TFile.Open(args.input)
    tree = f.Get("tree")
    max_ev = tree.GetEntries()
    if args.max_events > 0:
        max_ev = min(max_ev, args.max_events)

    print(f"Processing {max_ev} events from {sample}...")

    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("photon.7x11ClusterLr2E", 1)
    tree.SetBranchStatus("photon.7x11ClusterLr2Eta", 1)
    tree.SetBranchStatus("photon.weta2", 1)
    tree.SetBranchStatus("photon.unfudged_weta2", 1)
    tree.SetBranchStatus("photon_cluster.etaSample2", 1)
    tree.SetBranchStatus("photon_cluster.etamax2", 1)

    cell_vec = ROOT.std.vector('double')()
    eta_vec  = ROOT.std.vector('double')()
    tree.SetBranchAddress("photon.7x11ClusterLr2E",   cell_vec)
    tree.SetBranchAddress("photon.7x11ClusterLr2Eta",  eta_vec)

    h_fudged = ROOT.TH1F("h_weta2_fudged", ";w_{#eta_{2}};Events", 100, 0.005, 0.025)
    h_unfudged = ROOT.TH1F("h_weta2_unfudged", ";w_{#eta_{2}};Events", 100, 0.005, 0.025)
    h_computed = ROOT.TH1F("h_weta2_computed", ";w_{#eta_{2}};Events", 100, 0.005, 0.025)
    h_diff_comp_unfudged = ROOT.TH1F("h_weta2_diff_computed_unfudged",
                                      ";#Delta w_{#eta_{2}};Events", 100, -0.005, 0.005)
    h_diff_fudged_unfudged = ROOT.TH1F("h_weta2_diff_fudged_unfudged",
                                        ";#Delta w_{#eta_{2}};Events", 100, -0.005, 0.005)
    h_2d = ROOT.TH2F("h_weta2_2d", ";Unfudged w_{#eta_{2}};Computed w_{#eta_{2}}",
                      100, 0.005, 0.025, 100, 0.005, 0.025)

    arr_fudged = np.empty(max_ev, dtype=np.float64)
    arr_unfudged = np.empty(max_ev, dtype=np.float64)
    arr_computed = np.empty(max_ev, dtype=np.float64)
    n_good = 0

    leaf_fudged = tree.GetLeaf("photon.weta2")
    leaf_unfudged = tree.GetLeaf("photon.unfudged_weta2")
    leaf_etaSample2 = tree.GetLeaf("photon_cluster.etaSample2")
    leaf_etamax2 = tree.GetLeaf("photon_cluster.etamax2")

    n_skipped_coverage = 0
    for i in range(max_ev):
        tree.GetEntry(i)

        if cell_vec.size() != 77 or eta_vec.size() != 77:
            continue

        # Skip events with incomplete grid coverage (missing calorimeter cells)
        nz = sum(1 for j in range(77) if not (cell_vec[j] == 0.0 and eta_vec[j] == 0.0))
        if nz != 77:
            n_skipped_coverage += 1
            continue

        etaSample2 = leaf_etaSample2.GetValue()
        etamax2 = leaf_etamax2.GetValue()

        weta2_comp = ROOT.compute_weta2_corrected(cell_vec, eta_vec, etaSample2, etamax2)
        if weta2_comp < 0:
            continue

        wf = leaf_fudged.GetValue()
        wu = leaf_unfudged.GetValue()

        h_fudged.Fill(wf)
        h_unfudged.Fill(wu)
        h_computed.Fill(weta2_comp)
        h_diff_comp_unfudged.Fill(weta2_comp - wu)
        h_diff_fudged_unfudged.Fill(wf - wu)
        h_2d.Fill(wu, weta2_comp)

        arr_fudged[n_good] = wf
        arr_unfudged[n_good] = wu
        arr_computed[n_good] = weta2_comp
        n_good += 1

        if (i + 1) % 500000 == 0:
            print(f"  Processed {i+1}/{max_ev} events...")

    arr_fudged = arr_fudged[:n_good]
    arr_unfudged = arr_unfudged[:n_good]
    arr_computed = arr_computed[:n_good]

    print(f"\n{'='*60}")
    print(f"W_ETA_2 SUMMARY — {sample}")
    print(f"{'='*60}")
    print(f"Events analyzed: {n_good}")
    print(f"Events skipped (incomplete grid): {n_skipped_coverage}")

    diff_cu = arr_computed - arr_unfudged
    print(f"\nComputed vs Unfudged (all events):")
    print(f"  Mean diff:  {np.mean(diff_cu):.6f}")
    print(f"  Std diff:   {np.std(diff_cu):.6f}")
    print(f"  Max |diff|: {np.max(np.abs(diff_cu)):.6f}")

    # Outlier analysis: events where |diff| > 0.001 are caused by
    # CaloFillRectangularCluster missing cells that Athena's fresh
    # CaloCellList from the full CaloCellContainer captures.
    outlier_mask = np.abs(diff_cu) > 0.001
    n_outliers = np.sum(outlier_mask)
    print(f"\n  Outlier events (|diff| > 0.001): {n_outliers}/{n_good}")
    if n_outliers > 0:
        for idx in np.where(outlier_mask)[0]:
            print(f"    unfudged={arr_unfudged[idx]:.6f}  "
                  f"computed={arr_computed[idx]:.6f}  "
                  f"diff={diff_cu[idx]:.6f}")
        good = ~outlier_mask
        print(f"\n  Excluding outliers ({np.sum(good)} events):")
        print(f"    Mean diff:  {np.mean(diff_cu[good]):.6f}")
        print(f"    Std diff:   {np.std(diff_cu[good]):.6f}")
        print(f"    Max |diff|: {np.max(np.abs(diff_cu[good])):.6f}")

    diff_fu = arr_fudged - arr_unfudged
    print(f"\nFudged vs Unfudged:")
    print(f"  Mean diff:  {np.mean(diff_fu):.6f}")
    print(f"  Std diff:   {np.std(diff_fu):.6f}")

    fout = ROOT.TFile(args.output, "RECREATE")
    h_fudged.Write()
    h_unfudged.Write()
    h_computed.Write()
    h_diff_comp_unfudged.Write()
    h_diff_fudged_unfudged.Write()
    h_2d.Write()
    fout.Close()
    print(f"\nHistograms saved to {args.output}")


if __name__ == "__main__":
    main()
