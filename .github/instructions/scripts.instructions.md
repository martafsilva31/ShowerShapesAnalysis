---
applyTo: "scripts/**"
---

# Analysis Scripts

## Conventions

- Python scripts use `argparse` with `--input`, `--output`, `--max-events` flags.
- ROOT macros (`.C`) are run via `root -l -b -q 'script.C'`.
- Scripts read ntuples from `../ntuples/` (relative) or accept explicit `--input` paths.
- Output plots go to `../plots/`.

## Ntuple Structure

- Tree name: `tree`
- Key branches:
  - `photon.7x11ClusterLr2E` — `vector<double>`, 77 cells (7 eta × 11 phi, row-major)
  - `photon.reta`, `photon.unfudged_reta` — fudged / unfudged R_eta
  - `photon.rphi`, `photon.unfudged_rphi` — fudged / unfudged R_phi
  - `photon.weta2`, `photon.unfudged_weta2` — fudged / unfudged w_eta_2
  - `photon.pt`, `photon_cluster.eta2` — kinematics

## Shower Shape Definitions

- E277 = E(7η × 7φ) — central 7 phi cols (2–8) from all 7 eta rows
- E237 = E(3η × 7φ) — central 3 eta rows (2–4), central 7 phi cols (2–8)
- E233 = E(3η × 3φ) — central 3 eta rows (2–4), central 3 phi cols (4–6)
- R_eta = E237 / E277
- R_phi = E233 / E237
- w_eta_2 = sqrt( sum(E_i * eta_i^2)/sum(E_i) - (sum(E_i * eta_i)/sum(E_i))^2 )
  - Computed in a 3×5 (η×φ) window: eta rows 2–4, phi cols 3–7
  - Cell η spacing: 0.025 (Layer 2 granularity)
  - Uses relative positions (-1, 0, +1) × 0.025

## Cell-Energy Reweighting Method (ATL-COM-PHYS-2021-640, Section 5.2)

The reweighting method corrects both mean and variance of normalized cell energy distributions:
1. Normalized energy: `e_i = E_i / E_cluster`
2. Reweighted: `e_i^{RW} = (RMS_data_i / RMS_MC_i) × (e_i^MC − ē_MC_i) + ē_data_i`
3. Reweighted cell energy: `E_i^{RW} = (RMS_data_i / RMS_MC_i) × E_i^MC + (ē_data_i − (RMS_data/RMS_MC) × ē_MC_i) × E^MC`
4. Rescale to preserve cluster energy: `E_i^{RW} *= Σ E_i^MC / Σ E_i^{RW}`

Corrections are parameterized in bins of pT, |η|, conversion type, and optionally hit position (λ, quadrants).

Reference implementation in supervisor's code: `showershapereweighting/var.h` (`calcWeta`, `Energy` functions).

## Pipeline

The workflow is: compute script → histogram ROOT file (in `output/`) → plot macro (PDF in `output/plots/`).
- `compute_*.py` scripts read ntuples, compute shower shapes from cells, and write
  ROOT histogram files with fudged, unfudged, and cell-computed distributions.
- `plot_*.C` ROOT macros read those histogram files and produce PDF with
  an overlay panel (log-y) and a ratio-to-unfudged panel.
- Shared utilities in `utils.py` handle ntuple I/O, energy windows, and argparse.

## When Adding New Scripts

- Include docstring explaining the physics.
- Add `if __name__ == "__main__": main()` pattern.
- Update README task list.
