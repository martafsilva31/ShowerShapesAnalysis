#!/usr/bin/env python3
"""
Compute w_eta_2 from Layer 2 cell energies and compare with stored values.

Implements the same formula as var.h calcWeta:
    w_eta_2 = sqrt( sum(E_i * eta_i^2) / sum(E_i)
                  - (sum(E_i * |eta_i|) / sum(E_i))^2 )

where eta_i are approximate cell eta positions computed from
photon_cluster.eta2 + (ie - 1) * 0.025, and |eta_i| is used for the
first moment (matching the abs() in var.h calcWeta).

Uses a compiled C++ function for the inner computation to maximise speed.

Usage:
    python compute_weta_2.py -i ntuples/mc23e/mc23e_700770_Zeeg.root -o output/weta2_Zeeg.root
"""

import ROOT
import numpy as np
from utils import make_parser, extract_sample_label

ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>
// Matches var.h calcWeta:
//   sumEEta   += E * abs(eta_cell)
//   sumEEtaSq += E * eta_cell * eta_cell
// Cell etas approximated as: eta_cluster + (ie - 1) * 0.025
double compute_weta2_cpp(const std::vector<double>& cells, double eta_cluster) {
    if (cells.size() != 77) return -999.0;
    const int nphi = 11;
    const double cell_size = 0.025;
    // 3x5 window: eta rows [2,5), phi cols [3,8)
    double sumE = 0, sumEeta = 0, sumEeta2 = 0;
    for (int ie = 0; ie < 3; ie++) {
        double eta_cell = eta_cluster + (ie - 1) * cell_size;
        for (int ip = 3; ip < 8; ip++) {
            double e = cells[(ie + 2) * nphi + ip];
            sumE     += e;
            sumEeta  += e * std::abs(eta_cell);   // abs() as in var.h
            sumEeta2 += e * eta_cell * eta_cell;   // no abs, as in var.h
        }
    }
    if (sumE <= 0) return -999.0;
    double var = sumEeta2 / sumE - std::pow(sumEeta / sumE, 2);
    return (var >= 0) ? std::sqrt(var) : -999.0;
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
    tree.SetBranchStatus("photon.weta2", 1)
    tree.SetBranchStatus("photon.unfudged_weta2", 1)
    tree.SetBranchStatus("photon_cluster.eta2", 1)

    cell_vec = ROOT.std.vector('double')()
    tree.SetBranchAddress("photon.7x11ClusterLr2E", cell_vec)

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
    leaf_eta2 = tree.GetLeaf("photon_cluster.eta2")

    for i in range(max_ev):
        tree.GetEntry(i)

        if cell_vec.size() != 77:
            continue

        eta2 = leaf_eta2.GetValue()
        weta2_comp = ROOT.compute_weta2_cpp(cell_vec, eta2)
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

    diff_cu = arr_computed - arr_unfudged
    print(f"\nComputed vs Unfudged:")
    print(f"  Mean diff:  {np.mean(diff_cu):.6f}")
    print(f"  Std diff:   {np.std(diff_cu):.6f}")
    print(f"  Max |diff|: {np.max(np.abs(diff_cu)):.6f}")

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
