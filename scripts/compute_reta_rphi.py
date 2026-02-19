#!/usr/bin/env python3
"""
Compute Reta and Rphi from Layer 2 cell energies and compare with stored values.

Reta = E237 / E277 = E(3 eta × 7 phi) / E(7 eta × 7 phi) in Layer 2
Rphi = E233 / E237 = E(3 eta × 3 phi) / E(3 eta × 7 phi) in Layer 2

The 7×11 cluster in Layer 2 is stored as:
- 7 cells in eta direction × 11 cells in phi direction
- Flat vector of 77 cells in row-major order: [eta][phi]
- Cell index = eta_row * 11 + phi_col

The E277 uses the central 7 phi cells (columns 2-8) from all 7 eta rows.
The E237 uses the central 3 eta rows (rows 2-4) and central 7 phi columns (2-8).
The E233 uses the central 3 eta rows (rows 2-4) and central 3 phi columns (4-6).

NOTE: The computed Reta/Rphi matches the stored unfudged values when the photon
is well-centered in the cluster. Off-center photons or hot cells will show 
discrepancies because the stored values use the actual cluster center, not the 
grid center.
"""

import ROOT
import numpy as np
from glob import glob
import os

def get_energy_window(cell_energies, neta_cluster, nphi_cluster, neta_window, nphi_window):
    """
    Extract sum of energies in a centered window from a rectangular cluster.
    
    The 7x11 cluster in ATLAS Layer 2 is: 7 cells in eta × 11 cells in phi.
    Stored as flat array in row-major order: [eta][phi], so phi varies fastest.
    
    cell_energies: flat array of cell energies 
    neta_cluster, nphi_cluster: cluster dimensions (7, 11)
    neta_window, nphi_window: window dimensions to extract
    
    Returns the sum of energies in the centered window.
    """
    if len(cell_energies) != neta_cluster * nphi_cluster:
        return -999.0
    
    # Center indices (0-indexed)
    eta_c = neta_cluster // 2  # 3 for 7 cells
    phi_c = nphi_cluster // 2  # 5 for 11 cells
    
    # Calculate window boundaries (centered)
    eta_start = eta_c - neta_window // 2
    eta_end = eta_start + neta_window
    phi_start = phi_c - nphi_window // 2
    phi_end = phi_start + nphi_window
    
    # Sum over the window
    total = 0.0
    for ie in range(eta_start, eta_end):
        for ip in range(phi_start, phi_end):
            total += cell_energies[ie * nphi_cluster + ip]
    
    return total


def compute_shower_shapes(cell_energies_7x11):
    """
    Compute Reta and Rphi from 7x11 Layer 2 cell energies.
    
    Note: The 7x11 window in L2 means 7 cells in eta, 11 cells in phi.
    """
    neta, nphi = 7, 11
    
    # E(7x7) - use full eta range (7), central 7 in phi from 11
    E_7x7 = get_energy_window(cell_energies_7x11, neta, nphi, 7, 7)
    
    # E(3x7) - central 3 in eta, central 7 in phi
    E_3x7 = get_energy_window(cell_energies_7x11, neta, nphi, 3, 7)
    
    # E(3x3) - central 3 in eta, central 3 in phi
    E_3x3 = get_energy_window(cell_energies_7x11, neta, nphi, 3, 3)
    
    # Reta = E(3x7) / E(7x7)
    Reta = E_3x7 / E_7x7 if E_7x7 > 0 else -999.0
    
    # Rphi = E(3x3) / E(3x7)
    Rphi = E_3x3 / E_3x7 if E_3x7 > 0 else -999.0
    
    return Reta, Rphi, E_3x3, E_3x7, E_7x7


def analyze_file(filename, max_events=-1):
    """Analyze a single ROOT file."""
    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Error opening {filename}")
        return None
    
    tree = f.Get("tree")
    if not tree:
        print(f"No tree found in {filename}")
        return None
    
    results = {
        'reta_stored': [],
        'reta_unfudged': [],
        'reta_computed': [],
        'rphi_stored': [],
        'rphi_unfudged': [],
        'rphi_computed': [],
        'pt': [],
        'eta': [],
    }
    
    # Set up branch addresses for vector<double> branches
    cell_energies_vec = ROOT.std.vector('double')()
    tree.SetBranchAddress("photon.7x11ClusterLr2E", cell_energies_vec)
    
    nentries = tree.GetEntries()
    if max_events > 0:
        nentries = min(nentries, max_events)
    
    for i in range(nentries):
        tree.GetEntry(i)
        
        # Get stored values
        reta_stored = tree.GetLeaf("photon.reta").GetValue()
        reta_unfudged = tree.GetLeaf("photon.unfudged_reta").GetValue()
        rphi_stored = tree.GetLeaf("photon.rphi").GetValue()
        rphi_unfudged = tree.GetLeaf("photon.unfudged_rphi").GetValue()
        
        pt = tree.GetLeaf("photon.pt").GetValue()
        eta = tree.GetLeaf("photon_cluster.eta2").GetValue()
        
        # Get cell energies from 7x11 Layer 2 cluster (vector<double>)
        cell_energies = list(cell_energies_vec)
        
        if len(cell_energies) == 77:
            reta_comp, rphi_comp, E_3x3, E_3x7, E_7x7 = compute_shower_shapes(cell_energies)
            
            results['reta_stored'].append(reta_stored)
            results['reta_unfudged'].append(reta_unfudged)
            results['reta_computed'].append(reta_comp)
            results['rphi_stored'].append(rphi_stored)
            results['rphi_unfudged'].append(rphi_unfudged)
            results['rphi_computed'].append(rphi_comp)
            results['pt'].append(pt)
            results['eta'].append(eta)
    
    f.Close()
    return results


def compare_results(results):
    """Print comparison statistics."""
    reta_stored = np.array(results['reta_stored'])
    reta_unfudged = np.array(results['reta_unfudged'])
    reta_computed = np.array(results['reta_computed'])
    rphi_stored = np.array(results['rphi_stored'])
    rphi_unfudged = np.array(results['rphi_unfudged'])
    rphi_computed = np.array(results['rphi_computed'])
    
    print("\n" + "="*60)
    print("RETA COMPARISON")
    print("="*60)
    
    # Compare computed vs unfudged (should be very close)
    diff_unfudged = reta_computed - reta_unfudged
    print(f"\nComputed vs Unfudged (stored):")
    print(f"  Mean difference: {np.mean(diff_unfudged):.6f}")
    print(f"  Std difference:  {np.std(diff_unfudged):.6f}")
    print(f"  Max |diff|:      {np.max(np.abs(diff_unfudged)):.6f}")
    
    # Compare fudged vs unfudged
    diff_fudge = reta_stored - reta_unfudged
    print(f"\nFudged vs Unfudged (effect of fudging):")
    print(f"  Mean difference: {np.mean(diff_fudge):.6f}")
    print(f"  Std difference:  {np.std(diff_fudge):.6f}")
    
    print("\n" + "="*60)
    print("RPHI COMPARISON")
    print("="*60)
    
    diff_unfudged = rphi_computed - rphi_unfudged
    print(f"\nComputed vs Unfudged (stored):")
    print(f"  Mean difference: {np.mean(diff_unfudged):.6f}")
    print(f"  Std difference:  {np.std(diff_unfudged):.6f}")
    print(f"  Max |diff|:      {np.max(np.abs(diff_unfudged)):.6f}")
    
    diff_fudge = rphi_stored - rphi_unfudged
    print(f"\nFudged vs Unfudged (effect of fudging):")
    print(f"  Mean difference: {np.mean(diff_fudge):.6f}")
    print(f"  Std difference:  {np.std(diff_fudge):.6f}")
    
    print(f"\nTotal events analyzed: {len(reta_stored)}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compute Reta/Rphi from cell energies')
    parser.add_argument('--input', '-i', required=True, help='Input ROOT file or directory')
    parser.add_argument('--max-events', '-n', type=int, default=10000, help='Max events to process (-1 for all)')
    parser.add_argument('--output', '-o', help='Output ROOT file for histograms')
    args = parser.parse_args()
    
    # Collect input files
    if os.path.isdir(args.input):
        files = glob(os.path.join(args.input, '**/*.root'), recursive=True)
    else:
        files = [args.input]
    
    print(f"Processing {len(files)} file(s)...")
    
    all_results = {
        'reta_stored': [], 'reta_unfudged': [], 'reta_computed': [],
        'rphi_stored': [], 'rphi_unfudged': [], 'rphi_computed': [],
        'pt': [], 'eta': [],
    }
    
    events_per_file = args.max_events // len(files) if args.max_events > 0 else -1
    
    for i, f in enumerate(files[:10]):  # Limit to first 10 files for quick test
        print(f"  [{i+1}/{min(len(files),10)}] {os.path.basename(f)}")
        results = analyze_file(f, events_per_file)
        if results:
            for key in all_results:
                all_results[key].extend(results[key])
    
    compare_results(all_results)
    
    # Optionally save histograms
    if args.output:
        fout = ROOT.TFile(args.output, "RECREATE")
        
        h_reta_diff = ROOT.TH1F("h_reta_diff", "Computed - Unfudged Reta;#Delta R_{#eta};Events", 100, -0.01, 0.01)
        h_rphi_diff = ROOT.TH1F("h_rphi_diff", "Computed - Unfudged Rphi;#Delta R_{#phi};Events", 100, -0.01, 0.01)
        
        for reta_c, reta_u in zip(all_results['reta_computed'], all_results['reta_unfudged']):
            h_reta_diff.Fill(reta_c - reta_u)
        for rphi_c, rphi_u in zip(all_results['rphi_computed'], all_results['rphi_unfudged']):
            h_rphi_diff.Fill(rphi_c - rphi_u)
        
        h_reta_diff.Write()
        h_rphi_diff.Write()
        fout.Close()
        print(f"\nHistograms saved to {args.output}")


if __name__ == "__main__":
    main()
