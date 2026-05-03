# Data–MC Normalisation Mismatch Study

Investigates kinematic normalisation differences between data and MC for the Z→llγ sample.
Generates comparison plots and a LaTeX compendium.

## Files

| File | Purpose |
|------|---------|
| `plot_kinematics.C` | Kinematic variable comparisons (standard) |
| `plot_kinematics_abs.C` | Absolute-normalisation kinematic comparisons |
| `run_kinematics.sh` | Driver for data24 ntuples |
| `run_kinematics_phase2.sh` | Driver for phase-2 (updated) ntuples |
| `make_kinematics_compendium.py` | LaTeX compendium generator |

## Running

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
cd scripts/data_mc_normalisation_mismatches
bash run_kinematics.sh
```
