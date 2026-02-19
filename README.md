# ShowerShapesAnalysis

Analysis scripts and grid submission for the ATLAS photon shower shape reweighting study.

## Directory Structure

```
ShowerShapesAnalysis/
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
