#!/usr/bin/env python3
"""
Compute R_phi from Layer 2 cell energies and compare with stored values.

R_phi = E233 / E237 = E(3 eta x 3 phi) / E(3 eta x 7 phi) in Layer 2

Uses a compiled C++ function for the inner computation to maximise speed.

Usage:
    python compute_rphi.py -i ntuples/mc23e/mc23e_700770_Zeeg.root -o output/rphi_Zeeg.root
"""

import ROOT
import numpy as np
from utils import make_parser, extract_sample_label

ROOT.gInterpreter.Declare(r"""
#include <vector>
double compute_rphi_cpp(const std::vector<double>& cells) {
    if (cells.size() != 77) return -999.0;
    const int nphi = 11;
    // E(3eta x 7phi): eta rows [2,5), phi cols [2,9)
    double E_3x7 = 0;
    for (int ie = 2; ie < 5; ie++)
        for (int ip = 2; ip < 9; ip++)
            E_3x7 += cells[ie * nphi + ip];
    // E(3eta x 3phi): eta rows [2,5), phi cols [4,7)
    double E_3x3 = 0;
    for (int ie = 2; ie < 5; ie++)
        for (int ip = 4; ip < 7; ip++)
            E_3x3 += cells[ie * nphi + ip];
    return (E_3x7 > 0) ? E_3x3 / E_3x7 : -999.0;
}
""")


def main():
    parser = make_parser('Compute R_phi from cell energies')
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
    tree.SetBranchStatus("photon.rphi", 1)
    tree.SetBranchStatus("photon.unfudged_rphi", 1)

    cell_vec = ROOT.std.vector('double')()
    eta_vec  = ROOT.std.vector('double')()
    tree.SetBranchAddress("photon.7x11ClusterLr2E", cell_vec)
    tree.SetBranchAddress("photon.7x11ClusterLr2Eta", eta_vec)

    h_fudged = ROOT.TH1F("h_rphi_fudged", ";R_{#phi};Events", 100, 0.5, 1.05)
    h_unfudged = ROOT.TH1F("h_rphi_unfudged", ";R_{#phi};Events", 100, 0.5, 1.05)
    h_computed = ROOT.TH1F("h_rphi_computed", ";R_{#phi};Events", 100, 0.5, 1.05)
    h_diff_comp_unfudged = ROOT.TH1F("h_rphi_diff_computed_unfudged",
                                      ";#Delta R_{#phi};Events", 100, -0.02, 0.02)
    h_diff_fudged_unfudged = ROOT.TH1F("h_rphi_diff_fudged_unfudged",
                                        ";#Delta R_{#phi};Events", 100, -0.02, 0.02)
    h_2d = ROOT.TH2F("h_rphi_2d", ";Unfudged R_{#phi};Computed R_{#phi}",
                      100, 0.5, 1.05, 100, 0.5, 1.05)

    arr_fudged = np.empty(max_ev, dtype=np.float64)
    arr_unfudged = np.empty(max_ev, dtype=np.float64)
    arr_computed = np.empty(max_ev, dtype=np.float64)
    n_good = 0

    leaf_fudged = tree.GetLeaf("photon.rphi")
    leaf_unfudged = tree.GetLeaf("photon.unfudged_rphi")

    n_skipped_coverage = 0
    for i in range(max_ev):
        tree.GetEntry(i)

        if cell_vec.size() != 77:
            continue

        # Skip events with incomplete grid coverage (missing calorimeter cells)
        nz = sum(1 for j in range(77) if not (cell_vec[j] == 0.0 and eta_vec[j] == 0.0))
        if nz != 77:
            n_skipped_coverage += 1
            continue

        rphi_comp = ROOT.compute_rphi_cpp(cell_vec)
        if rphi_comp < 0:
            continue

        rf = leaf_fudged.GetValue()
        ru = leaf_unfudged.GetValue()

        h_fudged.Fill(rf)
        h_unfudged.Fill(ru)
        h_computed.Fill(rphi_comp)
        h_diff_comp_unfudged.Fill(rphi_comp - ru)
        h_diff_fudged_unfudged.Fill(rf - ru)
        h_2d.Fill(ru, rphi_comp)

        arr_fudged[n_good] = rf
        arr_unfudged[n_good] = ru
        arr_computed[n_good] = rphi_comp
        n_good += 1

        if (i + 1) % 500000 == 0:
            print(f"  Processed {i+1}/{max_ev} events...")

    arr_fudged = arr_fudged[:n_good]
    arr_unfudged = arr_unfudged[:n_good]
    arr_computed = arr_computed[:n_good]

    print(f"\n{'='*60}")
    print(f"R_PHI SUMMARY — {sample}")
    print(f"{'='*60}")
    print(f"Events analyzed: {n_good}")
    print(f"Events skipped (incomplete grid): {n_skipped_coverage}")

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
