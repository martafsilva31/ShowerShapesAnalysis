#!/usr/bin/env python3
"""
Event-by-event comparison of v5 (no reorderToGrid) vs v9 (with reorderToGrid).

Computes R_eta, R_phi, and w_eta_2 from Lr2 cell energies in both ntuples,
applies the crack cut |eta| not in (1.37, 1.52), and reports max differences.

Usage:
    python compare_v5_v9.py \
        --v5 ../ntuples/mc23e_egam3_v3/mc23e_egam3_v5.root \
        --v9 ../ntuples/mc23e_egam3_v3/mc23e_egam3_v9.root
"""

import ROOT
import numpy as np
import argparse

ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>

struct ShowerShapes {
    double reta;
    double rphi;
    double weta2;
};

ShowerShapes compute_ss(const std::vector<double>& cells) {
    ShowerShapes ss = {-999., -999., -999.};
    if (cells.size() != 77) return ss;
    const int nphi = 11;

    // R_eta = E(3x7) / E(7x7)
    double E_3x7 = 0;
    for (int ie = 2; ie < 5; ie++)
        for (int ip = 2; ip < 9; ip++)
            E_3x7 += cells[ie * nphi + ip];
    double E_7x7 = 0;
    for (int ie = 0; ie < 7; ie++)
        for (int ip = 2; ip < 9; ip++)
            E_7x7 += cells[ie * nphi + ip];
    ss.reta = (E_7x7 > 0) ? E_3x7 / E_7x7 : -999.;

    // R_phi = E(3x3) / E(3x7)
    double E_3x3 = 0;
    for (int ie = 2; ie < 5; ie++)
        for (int ip = 4; ip < 7; ip++)
            E_3x3 += cells[ie * nphi + ip];
    ss.rphi = (E_3x7 > 0) ? E_3x3 / E_3x7 : -999.;

    // w_eta_2 in 3x5 window
    double sumE = 0, sumEe = 0, sumEe2 = 0;
    for (int ie = 0; ie < 3; ie++) {
        double eta_rel = (ie - 1) * 0.025;
        for (int ip = 0; ip < 5; ip++) {
            double e = cells[(2 + ie) * nphi + (3 + ip)];
            sumE += e;
            sumEe += e * eta_rel;
            sumEe2 += e * eta_rel * eta_rel;
        }
    }
    if (sumE > 0) {
        double var = sumEe2 / sumE - (sumEe / sumE) * (sumEe / sumE);
        ss.weta2 = (var >= 0) ? std::sqrt(var) : -999.;
    }
    return ss;
}
""")


def main():
    parser = argparse.ArgumentParser(description='Compare v5 vs v9 shower shapes')
    parser.add_argument('--v5', required=True, help='v5 ntuple (no reorderToGrid)')
    parser.add_argument('--v9', required=True, help='v9 ntuple (with reorderToGrid)')
    args = parser.parse_args()

    f5 = ROOT.TFile.Open(args.v5)
    f9 = ROOT.TFile.Open(args.v9)
    t5 = f5.Get("tree")
    t9 = f9.Get("tree")

    n5 = t5.GetEntries()
    n9 = t9.GetEntries()
    print(f"v5: {n5} events, v9: {n9} events")
    assert n5 == n9, "Event count mismatch!"
    nev = n5

    # Setup branches for both trees
    for t in [t5, t9]:
        t.SetBranchStatus("*", 0)
        t.SetBranchStatus("photon.7x11ClusterLr2E", 1)
        t.SetBranchStatus("photon.7x11ClusterLr2Eta", 1)
        t.SetBranchStatus("photon.eta", 1)

    cells5 = ROOT.std.vector('double')()
    cells9 = ROOT.std.vector('double')()
    eta5_vec = ROOT.std.vector('double')()
    eta9_vec = ROOT.std.vector('double')()
    t5.SetBranchAddress("photon.7x11ClusterLr2E", cells5)
    t9.SetBranchAddress("photon.7x11ClusterLr2E", cells9)
    t5.SetBranchAddress("photon.7x11ClusterLr2Eta", eta5_vec)
    t9.SetBranchAddress("photon.7x11ClusterLr2Eta", eta9_vec)

    # Counters
    n_total = 0
    n_crack = 0
    n_no_crack = 0
    n_skip = 0

    # Max diffs (no crack cut)
    max_d_reta_all = 0.0
    max_d_rphi_all = 0.0
    max_d_weta2_all = 0.0

    # Max diffs (with crack cut)
    max_d_reta = 0.0
    max_d_rphi = 0.0
    max_d_weta2 = 0.0

    # Max diffs (crack only)
    max_d_reta_crack = 0.0
    max_d_rphi_crack = 0.0
    max_d_weta2_crack = 0.0

    # Cell-level
    n_cells_identical = 0
    n_cells_differ = 0

    for i in range(nev):
        t5.GetEntry(i)
        t9.GetEntry(i)

        if cells5.size() != 77 or cells9.size() != 77:
            n_skip += 1
            continue

        # Check coverage (skip zero-padded cells)
        nz5 = sum(1 for j in range(77) if not (cells5[j] == 0.0 and eta5_vec[j] == 0.0))
        nz9 = sum(1 for j in range(77) if not (cells9[j] == 0.0 and eta9_vec[j] == 0.0))
        if nz5 != 77 or nz9 != 77:
            n_skip += 1
            continue

        n_total += 1

        # Get photon eta
        eta = t5.GetLeaf("photon.eta").GetValue()
        abs_eta = abs(eta)
        is_crack = 1.37 < abs_eta < 1.52

        # Compare cells
        cells_match = all(cells5[j] == cells9[j] for j in range(77))
        if cells_match:
            n_cells_identical += 1
        else:
            n_cells_differ += 1

        # Compute shower shapes
        ss5 = ROOT.compute_ss(cells5)
        ss9 = ROOT.compute_ss(cells9)

        if ss5.reta < 0 or ss9.reta < 0:
            continue

        d_reta = abs(ss5.reta - ss9.reta)
        d_rphi = abs(ss5.rphi - ss9.rphi)
        d_weta2 = abs(ss5.weta2 - ss9.weta2)

        # All events
        max_d_reta_all = max(max_d_reta_all, d_reta)
        max_d_rphi_all = max(max_d_rphi_all, d_rphi)
        max_d_weta2_all = max(max_d_weta2_all, d_weta2)

        if is_crack:
            n_crack += 1
            max_d_reta_crack = max(max_d_reta_crack, d_reta)
            max_d_rphi_crack = max(max_d_rphi_crack, d_rphi)
            max_d_weta2_crack = max(max_d_weta2_crack, d_weta2)
        else:
            n_no_crack += 1
            max_d_reta = max(max_d_reta, d_reta)
            max_d_rphi = max(max_d_rphi, d_rphi)
            max_d_weta2 = max(max_d_weta2, d_weta2)

    print(f"\n{'='*65}")
    print(f"  v5 vs v9 SHOWER SHAPE COMPARISON")
    print(f"{'='*65}")
    print(f"Total events:        {nev}")
    print(f"Skipped (coverage):  {n_skip}")
    print(f"Analyzed:            {n_total}")
    print(f"  Non-crack:         {n_no_crack}")
    print(f"  Crack (1.37-1.52): {n_crack}")
    print(f"Cells identical:     {n_cells_identical}")
    print(f"Cells differ:        {n_cells_differ}")

    print(f"\n--- ALL events (no eta cut) ---")
    print(f"  Max |Δ R_eta|:   {max_d_reta_all:.15e}")
    print(f"  Max |Δ R_phi|:   {max_d_rphi_all:.15e}")
    print(f"  Max |Δ w_eta_2|: {max_d_weta2_all:.15e}")

    print(f"\n--- NON-CRACK events (|eta| < 1.37 or |eta| > 1.52) ---")
    print(f"  Max |Δ R_eta|:   {max_d_reta:.15e}")
    print(f"  Max |Δ R_phi|:   {max_d_rphi:.15e}")
    print(f"  Max |Δ w_eta_2|: {max_d_weta2:.15e}")

    print(f"\n--- CRACK events only (1.37 < |eta| < 1.52) ---")
    print(f"  Max |Δ R_eta|:   {max_d_reta_crack:.15e}")
    print(f"  Max |Δ R_phi|:   {max_d_rphi_crack:.15e}")
    print(f"  Max |Δ w_eta_2|: {max_d_weta2_crack:.15e}")

    if max_d_reta == 0.0 and max_d_rphi == 0.0 and max_d_weta2 == 0.0:
        print(f"\n✓ CONFIRMED: Non-crack shower shapes are IDENTICAL between v5 and v9.")
        print(f"  reorderToGrid() is redundant outside the crack region.")
    else:
        print(f"\n✗ NON-CRACK DIFFERENCES FOUND — investigate further.")


if __name__ == "__main__":
    main()
