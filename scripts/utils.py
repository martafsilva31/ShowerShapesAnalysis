#!/usr/bin/env python3
"""
Shared utilities for shower shape computation scripts.

Provides common I/O, cell energy window extraction, and CLI helpers
used by compute_reta.py, compute_rphi.py, and compute_weta_2.py.
"""

import ROOT
import argparse
import os
import re
import math


def make_parser(description):
    """Create a standard argparse parser with common flags."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input', '-i', required=True,
                        help='Input ROOT file (ntuple)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output ROOT file for histograms')
    parser.add_argument('--max-events', '-n', type=int, default=-1,
                        help='Max events to process (-1 for all)')
    return parser


def make_plot_parser(description):
    """Create a standard argparse parser for plot scripts."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input', '-i', required=True,
                        help='Input ROOT file with histograms')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file path (PDF/PNG)')
    return parser


def open_ntuple(filename):
    """
    Open a ROOT ntuple file and return (TFile, TTree).
    Raises RuntimeError if the file or tree cannot be opened.
    """
    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open file: {filename}")

    tree = f.Get("tree")
    if not tree:
        raise RuntimeError(f"No 'tree' found in {filename}")

    return f, tree


def setup_cell_branch(tree):
    """
    Set up the branch address for the 7x11 Layer 2 cell energies.
    Returns the std::vector<double> object bound to the branch.
    """
    cell_vec = ROOT.std.vector('double')()
    tree.SetBranchAddress("photon.7x11ClusterLr2E", cell_vec)
    return cell_vec


def get_scalar(tree, branch_name):
    """Read a scalar value from the current tree entry via GetLeaf."""
    leaf = tree.GetLeaf(branch_name)
    if not leaf:
        raise RuntimeError(f"Branch '{branch_name}' not found")
    return leaf.GetValue()


def get_energy_window(cell_energies, neta_cluster, nphi_cluster,
                      neta_window, nphi_window):
    """
    Extract sum of energies in a centered window from a rectangular cluster.

    The cluster is stored as a flat array in row-major order: [eta][phi],
    so phi varies fastest.

    Parameters
    ----------
    cell_energies : list or array
        Flat array of cell energies (length = neta_cluster * nphi_cluster).
    neta_cluster, nphi_cluster : int
        Full cluster dimensions (e.g. 7, 11).
    neta_window, nphi_window : int
        Dimensions of the centered window to extract.

    Returns
    -------
    float
        Sum of energies in the window, or -999 if input is invalid.
    """
    if len(cell_energies) != neta_cluster * nphi_cluster:
        return -999.0

    eta_c = neta_cluster // 2
    phi_c = nphi_cluster // 2

    eta_start = eta_c - neta_window // 2
    eta_end = eta_start + neta_window
    phi_start = phi_c - nphi_window // 2
    phi_end = phi_start + nphi_window

    total = 0.0
    for ie in range(eta_start, eta_end):
        for ip in range(phi_start, phi_end):
            total += cell_energies[ie * nphi_cluster + ip]

    return total


def compute_weta2_from_cells(cell_energies, neta_cluster=7, nphi_cluster=11,
                              eta_cell_size=0.025):
    """
    Compute w_eta_2 from 7x11 Layer 2 cell energies.

    w_eta_2 = sqrt( sum(E_i * eta_i^2) / sum(E_i)
                  - (sum(E_i * eta_i) / sum(E_i))^2 )

    Standard energy-weighted eta width (variance formula) in a 3x5
    (eta x phi) window centred on the geometric centre of the 7x11 cluster.

    Uses relative eta positions: (-1, 0, +1) * eta_cell_size.
    The variance is translation-invariant, so relative positions give
    the same result as absolute cell etas.

    All cells are included (negative energies from noise subtraction
    are not filtered out), matching Athena / supervisor code.

    Parameters
    ----------
    cell_energies : list
        Flat array of 77 cell energies (7 eta x 11 phi, eta-major order).
    neta_cluster, nphi_cluster : int
        Cluster dimensions.
    eta_cell_size : float
        Eta granularity of Layer 2 cells (default 0.025).

    Returns
    -------
    float
        w_eta_2 value, or -999 if computation fails.
    """
    if len(cell_energies) != neta_cluster * nphi_cluster:
        return -999.0

    # 3x5 window: 3 in eta, 5 in phi, centred on the cluster
    eta_c = neta_cluster // 2   # 3
    phi_c = nphi_cluster // 2   # 5

    neta_w, nphi_w = 3, 5
    eta_start = eta_c - neta_w // 2  # 2
    phi_start = phi_c - nphi_w // 2  # 3

    sum_e = 0.0
    sum_e_eta = 0.0
    sum_e_eta2 = 0.0

    for ie in range(neta_w):
        # Relative eta position: (-1, 0, +1) * cell_size
        eta_rel = (ie - neta_w // 2) * eta_cell_size
        for ip in range(nphi_w):
            idx = (eta_start + ie) * nphi_cluster + (phi_start + ip)
            e = cell_energies[idx]
            sum_e += e
            sum_e_eta += e * eta_rel
            sum_e_eta2 += e * eta_rel * eta_rel

    if sum_e <= 0:
        return -999.0

    variance = sum_e_eta2 / sum_e - (sum_e_eta / sum_e) ** 2
    if variance < 0:
        return -999.0

    return math.sqrt(variance)


def extract_sample_label(filename):
    """
    Extract a human-readable sample label from a filename.

    Examples:
        'mc23e_700770_Zeeg.root' → 'Zeeg'
        'mc23e_700771_Zmumug.root' → 'Zmumug'
    """
    base = os.path.basename(filename)
    # Try pattern: mc23e_DSID_NAME.root
    m = re.search(r'mc23e_\d+_(\w+)\.root', base)
    if m:
        return m.group(1)
    # Fallback: strip extension
    return os.path.splitext(base)[0]
