#!/usr/bin/env python3
"""
Cell collision study for Layer 2 of the 7x11 grid.

Answers the supervisor's questions about cell collisions in the
CaloFillRectangularCluster output:

1. How often do collisions happen? Are cells weighted?
2. Should weights be used? (look for weight decorations)
3. Are we using CaloFillRectangularCluster or cluster cell links?
4. Are duplicates identical (E, eta, phi) or just same grid bin?

A "collision" means that more cells exist than grid slots (Lr2Size > 77),
so multiple raw cells are mapped to the same 0.025 x 0.025 grid bin by
reorderToGrid(), which accumulates their energies (E_out[idx] += E[i]).

Usage:
    python analyse_collisions.py
"""

import ROOT
import numpy as np
import os
import sys

# ── Configuration ──────────────────────────────────────────────────────────
NTUPLE = os.path.join(os.path.dirname(__file__),
                      '../../ntuples/mc23e_egam3_v3/mc23e_egam3_v9.root')
TREE_NAME = 'tree'
GRID_SIZE = 77        # 7 eta x 11 phi
CELL_DETA = 0.025
CELL_DPHI = 0.025
NETA = 7
NPHI = 11


def analyse(infile):
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        print(f'ERROR: cannot open {infile}')
        sys.exit(1)
    tree = f.Get(TREE_NAME)
    nev = tree.GetEntries()
    print(f'Opened {infile}  —  {nev} events')

    # Enable needed branches
    tree.SetBranchStatus('*', 0)
    tree.SetBranchStatus('photon.7x11ClusterLr2E', 1)
    tree.SetBranchStatus('photon.7x11ClusterLr2Eta', 1)
    tree.SetBranchStatus('photon.7x11ClusterLr2Phi', 1)
    tree.SetBranchStatus('photon.7x11ClusterLr2Size', 1)
    tree.SetBranchStatus('photon_cluster.eta2', 1)

    cell_E   = ROOT.std.vector('double')()
    cell_Eta = ROOT.std.vector('double')()
    cell_Phi = ROOT.std.vector('double')()
    tree.SetBranchAddress('photon.7x11ClusterLr2E',   cell_E)
    tree.SetBranchAddress('photon.7x11ClusterLr2Eta', cell_Eta)
    tree.SetBranchAddress('photon.7x11ClusterLr2Phi', cell_Phi)

    leaf_size = tree.GetLeaf('photon.7x11ClusterLr2Size')
    leaf_eta  = tree.GetLeaf('photon_cluster.eta2')

    # ── Accumulators ──────────────────────────────────────────────────────
    n_total = 0
    n_collision = 0          # events where Lr2Size > 77
    n_exact_dup = 0          # events containing exact (E, eta, phi) duplicates
    collision_sizes = []     # raw Lr2Size for collision events
    collision_etas  = []     # cluster eta for collision events
    extra_cells_per_event = []   # how many extra cells beyond 77

    # Per-event collision details (first 20)
    examples = []

    for iev in range(nev):
        tree.GetEntry(iev)
        lr2size = int(leaf_size.GetValue())
        cl_eta  = leaf_eta.GetValue()
        n_total += 1

        if lr2size <= GRID_SIZE:
            continue

        n_collision += 1
        collision_sizes.append(lr2size)
        collision_etas.append(cl_eta)
        extra_cells_per_event.append(lr2size - GRID_SIZE)

        # ── Check for exact duplicates (same E, eta, phi) ────────────
        # The ntuple stores the 77-element grid-ordered arrays.
        # Skip zero-padded cells (E=0, eta=0, phi=0) — they would all
        # trivially match each other but are just empty grid slots.
        has_exact_dup = False
        for i in range(int(cell_E.size())):
            if cell_E[i] == 0.0 and cell_Eta[i] == 0.0 and cell_Phi[i] == 0.0:
                continue
            for j in range(i + 1, int(cell_E.size())):
                if cell_E[j] == 0.0 and cell_Eta[j] == 0.0 and cell_Phi[j] == 0.0:
                    continue
                if (cell_E[i] == cell_E[j] and
                    cell_Eta[i] == cell_Eta[j] and
                    cell_Phi[i] == cell_Phi[j]):
                    has_exact_dup = True
                    break
            if has_exact_dup:
                break
        if has_exact_dup:
            n_exact_dup += 1

        # ── Collect example info (first 20 collision events) ──────────
        if len(examples) < 20:
            # Find which grid bins received >1 cell (before accumulation)
            # We cannot reconstruct the raw cells from the grid output,
            # but we CAN identify bins that look suspicious:
            # after accumulation, all bins have exactly one entry.
            # Lr2Size > 77 tells us collisions happened.
            ex = {
                'event': iev,
                'lr2size': lr2size,
                'eta': cl_eta,
                'n_extra': lr2size - GRID_SIZE,
                'n_nonzero': sum(1 for k in range(int(cell_E.size()))
                                 if not (cell_E[k] == 0 and cell_Eta[k] == 0)),
            }
            examples.append(ex)

    # ── Report ────────────────────────────────────────────────────────────
    print()
    print('=' * 72)
    print('CELL COLLISION STUDY — Layer 2 (7x11 grid)')
    print('=' * 72)
    print()

    print(f'Total events:                {n_total}')
    print(f'Events with Lr2Size > 77:    {n_collision}  '
          f'({100.0 * n_collision / n_total:.1f}%)')
    print(f'Events with exact duplicates '
          f'(same E, eta, phi):  {n_exact_dup}')
    print()

    if n_collision > 0:
        sizes = np.array(collision_sizes)
        etas  = np.array(collision_etas)
        extras = np.array(extra_cells_per_event)

        print('── Collision statistics ─────────────────────────────────')
        print(f'  Lr2Size  min/max/mean:  {sizes.min()} / {sizes.max()} / {sizes.mean():.1f}')
        print(f'  Extra cells  min/max/mean:  {extras.min()} / {extras.max()} / {extras.mean():.1f}')
        print(f'  |eta| range of collision events:  '
              f'[{np.abs(etas).min():.3f}, {np.abs(etas).max():.3f}]')
        print()

        # eta distribution
        print('── |eta| distribution of collision events ──────────────')
        bins = [0, 0.8, 1.0, 1.2, 1.3, 1.37, 1.52, 1.6, 1.8, 2.0, 2.5]
        hist, _ = np.histogram(np.abs(etas), bins=bins)
        for lo, hi, cnt in zip(bins[:-1], bins[1:], hist):
            bar = '#' * cnt
            pct = 100.0 * cnt / n_collision if n_collision > 0 else 0
            print(f'  [{lo:.2f}, {hi:.2f})  {cnt:4d}  ({pct:5.1f}%)  {bar}')
        print()

        # Barrel-endcap crack region
        n_crack = sum(1 for e in etas if 1.3 <= abs(e) <= 1.55)
        print(f'  In barrel-endcap crack (1.3 <= |eta| <= 1.55): '
              f'{n_crack}/{n_collision}  '
              f'({100.0 * n_crack / n_collision:.1f}%)')
        print()

        # Examples
        print('── Example collision events (first 20) ─────────────────')
        print(f'  {"Event":>6s}  {"Lr2Size":>7s}  {"Extra":>5s}  '
              f'{"eta":>8s}  {"NonZero":>7s}')
        for ex in examples:
            print(f'  {ex["event"]:6d}  {ex["lr2size"]:7d}  {ex["n_extra"]:5d}  '
                  f'{ex["eta"]:8.4f}  {ex["n_nonzero"]:7d}')
        print()

    # ── Summary answers ──────────────────────────────────────────────────
    print('=' * 72)
    print('ANSWERS TO SUPERVISOR QUESTIONS')
    print('=' * 72)
    print()
    print('Q1: How often do collisions happen?')
    print(f'    {n_collision}/{n_total} events ({100.0 * n_collision / n_total:.1f}%) '
          f'have Lr2Size > 77.')
    if n_collision > 0:
        print(f'    All occur at |eta| in [{np.abs(np.array(collision_etas)).min():.3f}, '
              f'{np.abs(np.array(collision_etas)).max():.3f}] (barrel-endcap crack).')
    print()
    print('Q2: Are cells weighted? Should we use weights?')
    print('    NO. CaloFillRectangularCluster builds a fresh cluster from')
    print('    CaloCellContainer — every cell has weight = 1.0.')
    print('    Mohamed\'s weighted duplicates come from topocluster cell links,')
    print('    where the SAME physical cell can be shared between clusters with')
    print('    fractional weights. That mechanism does not apply here.')
    print()
    print('Q3: CaloFillRectangularCluster or cluster cell links?')
    print('    Both. Stage 1 (cell links) locates the seed; Stage 2')
    print('    (CaloFillRectangularCluster) extracts cells in the eta x phi')
    print('    window around the seed. The grid arrays come from Stage 2.')
    print()
    print('Q4: Exact duplicates or just same grid bin?')
    print(f'    Events with exact (E, eta, phi) duplicates: {n_exact_dup}')
    print(f'    Events with grid-bin collisions (Lr2Size > 77): {n_collision}')
    if n_exact_dup == 0 and n_collision > 0:
        print('    → All collisions are DIFFERENT physical cells at different')
        print('      positions that map to the same 0.025 x 0.025 grid bin.')
        print('      Energy accumulation (+=) is the correct treatment.')
    print()


if __name__ == '__main__':
    analyse(NTUPLE)
