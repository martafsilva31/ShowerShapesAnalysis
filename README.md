# ShowerShapesAnalysis

Analysis scripts and grid submission for the ATLAS photon shower shape reweighting study.

## Directory Structure

```
ShowerShapesAnalysis/
├── .github/            # Development conventions
├── scripts/            # Analysis scripts
│   ├── utils.py               # Shared I/O and cell-energy utilities
│   ├── compute_reta.py        # Compute R_eta from 7×11 cells
│   ├── compute_rphi.py        # Compute R_phi from 7×11 cells
│   ├── compute_weta_2.py      # Compute w_eta_2 from 7×11 cells
│   ├── plot_reta.C            # ROOT macro: R_eta overlay + ratio
│   ├── plot_rphi.C            # ROOT macro: R_phi overlay + ratio
│   ├── plot_weta_2.C          # ROOT macro: w_eta_2 overlay + ratio
│   ├── closure_test/          # MC closure test pipeline
│   │   ├── config.h           # Shared config: branch names, geometry, helpers
│   │   ├── create_pseudodata.C
│   │   ├── derive_corrections.C
│   │   ├── apply_corrections.C
│   │   ├── validate_closure.C
│   │   ├── plot_closure.C
│   │   ├── run_closure_test.sh
│   │   └── run_closure_suite.sh
│   └── data_mc/               # Data-MC cell-energy reweighting (current)
│       ├── config.h                # Shared config: cuts, branches, geometry, formulas
│       ├── fill_histograms.C       # Two-pass pipeline: accumulate → correct → fill
│       ├── plot_shower_shapes.C    # Shower shape plots: validation, fudge, comparison, per-eta
│       ├── plot_cell_profiles.C    # 7×11 cell heatmaps + correction vector plots
│       ├── extract_chi2.C          # Extract chi-squared tables → report/chi2_*.tex
│       ├── run.sh                  # Driver: compile, run fill + plot scripts
│       └── old/                    # Archived previous-iteration scripts
├── grid/               # Grid submission scripts (pathena)
│   ├── samples/        # Sample lists (dataset names)
│   ├── download_ntuples.sh    # Download + merge ntuples from grid
│   ├── submit_data22.sh
│   ├── submit_data24.sh
│   └── submit_mc23e_Zllg.sh
├── configs/            # Configuration files
├── ntuples/            # Output ntuples (gitignored)
│   ├── data22/         # Data 2022 (all 0 events — GRL mismatch, see below)
│   ├── data24/         # Data 2024
│   └── mc23e/          # MC23e (Zeeg 13 GB, Zmumug 23 GB)
├── output/             # Computed histograms + plots (gitignored)
│   ├── old/            # Archived output from previous iteration
│   └── cell_energy_reweighting_Francisco_method/
│       └── data24/{channel}/{scenario}/
│           ├── histograms.root
│           └── plots/  # 30 PDF comparison plots
├── report/             # Reports and documentation
│   ├── egam3_problem_report.md       # DAOD_EGAM3 problem analysis
│   ├── weta2_investigation_summary.md # w_eta_2 study summary
│   ├── closure_test_report.tex        # MC closure test report
│   ├── cell_reweighting_report.tex    # Physics report: M1/M2 pipeline (current)
│   ├── chi2_yields.tex               # Auto-generated event yield table
│   ├── chi2_summary.tex              # Auto-generated chi-squared summary
│   └── chi2_tables.tex               # Auto-generated per-eta chi-squared tables
└── setup.sh            # Environment setup
```

## Quick Start

```bash
# 1. Source the environment (from the ShowerShapes/ directory)
cd /project/atlas/users/mfernand/QT/ShowerShapes
source ShowerShapesAnalysis/setup.sh

# 2. Compute shower shapes (one script per variable, one run per sample)
cd ShowerShapesAnalysis/scripts
python compute_reta.py   -i ../ntuples/mc23e/mc23e_700770_Zeeg.root -o ../output/reta_Zeeg.root
python compute_rphi.py   -i ../ntuples/mc23e/mc23e_700770_Zeeg.root -o ../output/rphi_Zeeg.root
python compute_weta_2.py -i ../ntuples/mc23e/mc23e_700770_Zeeg.root -o ../output/weta2_Zeeg.root

# 3. Generate plots (ROOT macros — reads histogram files, outputs PDF)
root -l -b -q 'plot_reta.C("../output/reta_Zeeg.root","../output/plots/reta_Zeeg")'
root -l -b -q 'plot_rphi.C("../output/rphi_Zeeg.root","../output/plots/rphi_Zeeg")'
root -l -b -q 'plot_weta_2.C("../output/weta2_Zeeg.root","../output/plots/weta2_Zeeg")'
```

## MC Closure Test Pipeline

Validates cell-energy reweighting corrections using synthetic distortions:

```bash
cd scripts/closure_test
./run_closure_test.sh        # Single distortion level
./run_closure_suite.sh       # Suite: 1%, 3%, 5% distortions
```

See `scripts/closure_test/README.md` for full documentation.

## Data-MC Comparison Pipeline (Legacy)

Previous-iteration scripts that compared data and MC shower shapes at three
selection levels. These scripts have been archived to `scripts/data_mc/old/`
and superseded by the cell-energy reweighting pipeline above.

## Cell-Energy Reweighting Pipeline (Data vs MC)

Derives and applies per-cell energy corrections to photon shower shapes
(R_eta, R_phi, w_eta_2) using two methods:
- **M1 (flat shift)**: $E'_k = E_k + \Delta_k \times E_\mathrm{total}$
- **M2 (shift+stretch)**: $E'_k = E_\mathrm{total} \times \mathrm{shift}_k + \mathrm{stretch}_k \times E_k$

M2 is equivalent to Francisco's `photoncellbasedrw` method.

**Channels**: `eegamma` (Z→eeγ), `mumugamma` (Z→μμγ), `llgamma` (combined)
**Conversion scenarios**: `baseline` (unconverted), `converted`, `all_conv` (inclusive)

Selection: pT > 10 GeV, |η| < 2.37, crack excluded, loose isolation,
mll ∈ [40, 83] GeV, mllg ∈ [80, 100] GeV, ΔR(lep, γ) > 0.4.
MC weight: w = w_MC × w_μ × σ.

**Architecture**: Two-pass C++ pipeline compiled as standalone executable.
Pass 1 accumulates data/MC cell statistics, computes corrections.
Pass 2 applies M1 and M2 corrections and fills comparison histograms.

```bash
# Set up ROOT
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
asetup Athena,25.0.40 --quiet

# Run single channel/scenario
cd scripts/data_mc
./run.sh eegamma baseline

# Run all 9 combinations (3 channels × 3 scenarios)
./run.sh --batch

# Plot-only (reuses existing histograms.root)
./run.sh eegamma baseline --plot-only

# Extract chi-squared tables for the report
root -l -b -q 'extract_chi2.C'
```

Output: `output/cell_energy_reweighting_Francisco_method/data24/{channel}/{scenario}/`
containing `histograms.root` and `plots/` directory.

**Scenario matrix** (all 9 combinations):

| Channel | baseline | converted | all_conv |
|---------|:--------:|:---------:|:--------:|
| eegamma | ✅ 78K data / 956K MC | ✅ 16K / 240K | ✅ 94K / 1.2M |
| mumugamma | ✅ 127K / 1.8M | ✅ 29K / 516K | ✅ 156K / 2.4M |
| llgamma | ✅ 205K / 2.8M | ✅ 46K / 756K | ✅ 251K / 3.5M |

PDFs produced per scenario:
- `rew_{reta,rphi,weta2}.pdf` — SET B: cell-computed Data vs MC, M1, M2 (13 pages/eta)
- `computed_vs_stored.pdf` — cell-computed vs branch validation
- `fudge_factors.pdf` — SET A: Data vs MC unfudged vs fudged
- `cell_{data,mc,mc_m1,mc_m2}.pdf` — 7×11 cell fraction heatmaps
- `cell_shift.pdf`, `cell_stretch.pdf` — M2 correction vectors
- `cell_profiles_*.pdf` — before/after heatmap comparisons (baseline only)

## Grid Submission

All grid scripts include full environment setup (setupATLAS, asetup, lsetup, voms-proxy).
Submit from the `NTupleMaker_workspace/run/` directory:

```bash
bash ShowerShapesAnalysis/grid/submit_data24.sh      # Data 2024 run 473235
bash ShowerShapesAnalysis/grid/submit_mc23e_Zllg.sh   # MC23e Z→llγ
```

## Related Repositories

| Repo | Location | Remote |
|------|----------|--------|
| **NTupleMaker** | `../NTupleMaker_workspace/source/NTupleMaker/` | `gitlab.cern.ch/femarta/cellntuplemaker` |
| **showershapereweighting** | `../showershapereweighting/` | `gitlab.cern.ch/mobelfki/showershapereweighting` |

## Known Issues

- **Data22 has 0 events**: Runs 428648–429027 (periods A–E) are not covered by the 2022 GRL, which starts at run 431810 (period F). Need to re-download runs from periods F/H/J.

## Developer Workflow

### Adding a New Script

1. Place it in `scripts/`.
2. Include a docstring or header comment describing inputs/outputs.
3. Accept input file paths as CLI arguments (no hardcoded paths).
4. Follow the naming convention: `<verb>_<observable>.py` (e.g., `compute_weta_2.py`).
5. Use shared utilities from `utils.py` for file I/O and cell-energy extraction.

### Adding a New Grid Dataset

1. Verify the run is in the GRL (see `.github/instructions/grid.instructions.md`).
2. Create `grid/submit_<dataset>.sh` following the existing template.
3. Add sample dataset names to `grid/samples/`.


