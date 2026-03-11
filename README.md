# ShowerShapesAnalysis

Analysis scripts and grid submission for the ATLAS photon shower shape reweighting study.

## Directory Structure

```
ShowerShapesAnalysis/
├── .github/            # Copilot instructions (global + scoped)
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
│   └── data_mc/               # Data-MC comparison + reweighting
│       ├── compute_and_compare.C   # Initial data-vs-MC comparison
│       ├── validate_data_mc.C      # Validation (4 samples, chi²)
│       ├── plot_data_mc.C          # Final overlay plots
│       └── run_data_mc.sh          # Driver script (6-step pipeline)
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
│   ├── *.root          # Histogram ROOT files from compute scripts
│   ├── plots/          # PDF plots from ROOT macros
│   └── data_mc_comparison/  # Data-MC pipeline output
├── report/             # Reports and documentation
│   ├── egam3_problem_report.md       # DAOD_EGAM3 problem analysis
│   ├── weta2_investigation_summary.md # w_eta_2 study summary
│   └── data_mc_comparison_report.tex  # Data-MC comparison report
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

## Data-MC Comparison Pipeline

```bash
# Download ntuples from completed grid jobs
cd grid && ./download_ntuples.sh && cd ..

# Run the full pipeline (compute, derive corrections, apply, validate, plot)
cd scripts/data_mc && ./run_data_mc.sh
```

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

### Copilot / AI-Assisted Development

This repo includes `.github/copilot-instructions.md` and scoped instruction files in
`.github/instructions/` that give Copilot context about:
- The physics (shower shapes, EFT, cell-level ntuples)
- The software stack (Athena, ROOT, pathena)
- Conventions for scripts and grid jobs

These files are automatically picked up by GitHub Copilot in VS Code.

## Task List

- [x] Repository structure and cleanup
- [x] Merge split ntuples (data22, mc23e)
- [x] Diagnose data22 GRL issue
- [x] Prepare data24 grid submission
- [x] GitHub remote setup
- [x] Copilot instructions
- [x] Add fudged variable on top of R_eta / R_phi plots
- [x] w_eta_2 comparison plot (compute + fudged/unfudged overlay)
- [x] MC closure test pipeline (cell-energy reweighting)
- [x] Data-MC comparison pipeline
- [x] Add L2 cell eta/phi branches to NTupleMaker
- [x] Grid production v3 (task 49065711, 35M events, AODs with L2 eta/phi branches)
- [ ] Download v3 grid output and run analysis
- [ ] Resolve w_eta_2 computation discrepancy (see report/weta2_investigation_summary.md)
- [ ] Apply Francisco's shower shape reweighting correction
