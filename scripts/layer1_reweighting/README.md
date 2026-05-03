# Layer 1 Cell-Energy Reweighting Pipeline

Reweights EM calorimeter Layer 1 cell energy fractions in MC to match data, applying M1 and M2
corrections in bins of photon η (and optionally pT).

## Files

| File | Purpose |
|------|---------|
| `config.h` | Binning constants |
| `fill_histograms.C` | Two-pass fill: accumulate fractions → compute corrections → apply to MC |
| `plot_shower_shapes.C` | Shower shape comparison plots |
| `plot_cell_profiles.C` | Cell heatmaps and correction maps |
| `extract_chi2.C` | Chi-squared tables |
| `run_layer1_final.sh` | Full pipeline driver |

## Running

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
cd scripts/layer1_reweighting
bash run_layer1_final.sh
```
