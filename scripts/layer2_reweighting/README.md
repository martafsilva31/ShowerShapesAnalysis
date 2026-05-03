# Layer 2 Cell-Energy Reweighting Pipeline

Reweights EM calorimeter Layer 2 cell energy fractions in MC to match data, in bins of photon η
(and optionally pT or pileup μ). Computes M1 (shift/stretch) and M2 (delta) corrections per cell.

## Files

| File | Purpose |
|------|---------|
| `config.h` | Binning constants (η, pT, cluster size) |
| `fill_histograms.C` | Two-pass fill: accumulate fractions → compute corrections → apply to MC |
| `plot_shower_shapes.C` | Shower shape comparison plots (data vs corrected/uncorrected MC) |
| `plot_cell_profiles.C` | 7×11 cell heatmaps and correction maps |
| `plot_mu_scaling.C` | Pileup (μ) scaling diagnostic plots |
| `extract_chi2.C` | Chi-squared tables for shower shape compatibility |
| `extract_comparison_final.C` | Side-by-side variant comparison |
| `run_layer2_final.sh` | Full 6-phase pipeline driver |
| `make_compendiums.py` | LaTeX compendium generator for all variants/scenarios |

## Variants

| Variant | Binning |
|---------|---------|
| `eta_loose` | η only, loose isolation |
| `eta_tight` | η only, tight isolation |
| `eta_pt_loose` | η × pT, loose isolation |
| `eta_pt_tight` | η × pT, tight isolation |

## Running

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
cd scripts/layer2_reweighting
bash run_layer2_final.sh            # all variants, all scenarios
bash run_layer2_final.sh eta_loose  # single variant
```

## Output

Results written to `output/Layer_2/{variant}/llgamma/{scenario}/` where scenario ∈ {inclusive, converted, unconverted}.

Final compendium PDFs (tracked in git) are at `output/Layer_2/{variant}/llgamma/{scenario}/result_compendium_{variant}_llgamma_{scenario}.pdf`.
