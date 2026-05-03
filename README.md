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
│       ├── config.h                # Shared config: cuts, branches, geometry, formulas, pT/mu bins
│       ├── fill_histograms.C       # Two-pass pipeline: accumulate → correct → fill (eta/eta_pt/eta_mu × loose/tight iso)
│       ├── plot_shower_shapes.C    # Shower shape plots: per-eta PDFs + per-pT/per-mu PDFs
│       ├── plot_cell_profiles.C    # 7×11 cell heatmaps + correction vector plots
│       ├── plot_mu_scaling.C       # Pileup study: χ²/ndf and M1 shift vs ⟨μ⟩ (eta_mu variants)
│       ├── extract_chi2.C          # Extract chi-squared tables for a single variant
│       ├── extract_comparison_final.C  # Cross-variant chi2 comparison → report/chi2_variant_comparison.tex
│       ├── run_layer2_final.sh     # Full pipeline driver: 5 variants (eta/eta_pt/eta_mu × loose/tight)
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
│   └── Layer_2/        # Current output: 5 variants × 3 scenarios
│       ├── make_compendiums.py    # Generates LaTeX compendium PDFs for all variants
│       ├── eta_loose/             # η-only binning, loose isolation
│       ├── eta_tight/             # η-only binning, tight isolation
│       ├── eta_pt_loose/          # η×pT binning (14×6 bins), loose isolation
│       ├── eta_pt_tight/          # η×pT binning (14×6 bins), tight isolation
│       └── eta_mu_loose/          # η×⟨μ⟩ binning (14×4 bins), loose isolation — pileup study
│           └── {channel}/{scenario}/
│               ├── histograms.root
│               └── plots/  # shower shape PDFs, cell heatmaps, fudge factor plots
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

## Cell-Energy Reweighting Pipeline (Layer 2, Data vs MC)

Derives and applies per-cell energy corrections to photon shower shapes
(R_eta, R_phi, w_eta_2) using two methods:
- **M1 (flat shift)**: $E'_k = E_k + \Delta_k \times E_\mathrm{total}$
- **M2 (shift+stretch)**: $E'_k = E_\mathrm{total} \times \mathrm{shift}_k + \mathrm{stretch}_k \times E_k$

M2 is equivalent to Francisco's `photoncellbasedrw` method and is the **recommended** correction.

**Channel**: `llgamma` (Z→eeγ + Z→μμγ combined)
**Conversion scenarios**: `unconverted`, `converted`, `inclusive` (hadd of unc+conv)

Five **pipeline variants** are studied:

| Variant | Binning | Isolation |
|---------|---------|--------|
| `eta_loose` | 14 η bins | Loose |
| `eta_tight` | 14 η bins | Tight |
| `eta_pt_loose` | 14 η × 6 pT bins | Loose (**recommended**) |
| `eta_pt_tight` | 14 η × 6 pT bins | Tight |
| `eta_mu_loose` | 14 η × 4 ⟨μ⟩ bins | Loose (pileup study) |

pT bins (GeV): [10, 15, 20, 25, 30, 40, 1000].
⟨μ⟩ bins: [0, 40, 55, 70, 120].

Selection: pT > 10 GeV, |η| < 2.37, crack [1.37,1.52] excluded, loose/tight isolation,
mll ∈ [40, 83] GeV, mllg ∈ [80, 100] GeV, ΔR(lep, γ) > 0.4. No photon ID.

**Architecture**: Two-pass C++ pipeline run via ROOT.
Pass 1 accumulates data/MC cell statistics, computes M1/M2 corrections.
Pass 2 applies corrections and fills comparison histograms.

```bash
# Set up ROOT (LCG)
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

# Run all 4 variants (fill + plot + chi2 + compendium for all 3 scenarios each)
cd scripts/data_mc
bash run_layer2_final.sh

# Or run a single variant
VARIANTS="eta_loose" bash run_layer2_final.sh

# Extract cross-variant chi2 comparison table → report/chi2_variant_comparison.tex
root -l -b -q 'extract_comparison_final.C()'

# Regenerate all compendium PDFs
cd ../../output/Layer_2
python3 make_compendiums.py
```

Output per variant/scenario: `output/Layer_2/{variant}/{channel}/{scenario}/`
containing `histograms.root` and `plots/`.

Compendium PDFs: `output/Layer_2/{variant}/{channel}/{scenario}/result_compendium_{variant}_{channel}_{scenario}.pdf`

PDFs produced per scenario:
- `rew_{reta,rphi,weta2}.pdf` — per-eta: Data vs MC, M1, M2 (14 pages/eta bin)
- `rew_integrated.pdf` — all-eta integrated (3 pages)
- `rew_{var}_pt{PP}.pdf` — per-pT shower shapes (eta_pt variants only, 6 PDFs × 3 vars)
- `rew_{var}_mu{PP}.pdf` — per-⟨μ⟩ shower shapes (eta_mu variants only, 4 PDFs × 3 vars)
- `mu_chi2_scaling.pdf` — χ²/ndf vs ⟨μ⟩ bin (eta_mu variants only)
- `mu_m1_shift_scaling.pdf` — mean |M1 shift| vs ⟨μ⟩ bin (eta_mu variants only)
- `computed_vs_stored.pdf`, `computed_vs_stored_eta.pdf` — cell-computed vs branch validation
- `fudge_factors.pdf`, `fudge_factors_eta.pdf` — fudge factor comparison
- `cell_{data,mc,mc_m1,mc_m2}.pdf` — 7×11 cell fraction heatmaps
- `cell_shift.pdf`, `cell_stretch.pdf` — M2 correction vectors

**Event yields** (llgamma, loose isolation): ~205k data / ~2.8M MC (unconverted), ~46k / ~756k (converted).

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


