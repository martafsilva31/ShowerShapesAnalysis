# Shower Shape Analysis

Analysis scripts and grid submission for the ATLAS photon shower shape reweighting study.

## Directory Structure

```
analysis/
├── scripts/        # Analysis scripts (compute_reta_rphi.py, plot_reta_rphi.C, ...)
├── grid/           # Grid submission scripts (pathena)
│   └── samples/    # Sample lists (dataset names)
├── configs/        # Configuration files
├── ntuples/        # Output ntuples (gitignored)
│   ├── data22/     # Data 2022 ntuples
│   ├── data24/     # Data 2024 ntuples (future)
│   └── mc23e/      # MC23e ntuples (Zeeg, Zmumug)
├── plots/          # Output plots (gitignored)
└── setup.sh        # Environment setup
```

## Setup

```bash
source setup.sh
```

## Related Repositories

- **NTupleMaker**: `https://gitlab.cern.ch/femarta/cellntuplemaker.git`
  - Located at `../NTupleMaker_workspace/source/NTupleMaker/`
- **Shower Shape Reweighting**: `https://gitlab.cern.ch/mobelfki/showershapereweighting.git`
  - Located at `../showershapereweighting/`

## Task List

- [x] Compute r_eta and r_phi from Layer 2 cells (7x7 window) and validate against unfudged values
- [ ] Add fudged variable overlay on plots
- [ ] Make comparison plot for w_eta_2
- [ ] Understand why 2022 data did not pass the good run list
- [ ] Run on data24 (`data24_13p6TeV:data24_13p6TeV.00473235.physics_Main.merge.AOD.r15810_p6304`)
- [ ] Add random noise to MC cells
- [ ] Apply Francisco's cell-based shower shape reweighting correction
