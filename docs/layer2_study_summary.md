# Layer 2 Cell-Energy Reweighting — Study Summary

> **Status**: Stable release (April 2026, commit `bd6b1f4`).
> This document summarises the completed Layer 2 analysis and serves as context for the Layer 1 extension.

---

## 1. Physics Goal

Derive data-driven corrections for photon shower shapes ($R_\eta$, $R_\phi$, $w_{\eta 2}$) by comparing data and MC at the **cell level** in the second electromagnetic calorimeter layer (EM Layer 2). The cell-level approach (rather than shape-by-shape fudge factors) simultaneously improves all shower shapes built from those cells.

Reference: ATL-COM-PHYS-2021-640, Section 5.

---

## 2. Dataset

| Sample | Location | Notes |
|--------|----------|-------|
| Data 2024 | `/dcache/atlas/mfernand/qt_ntuples/data24/egam3.root` (Z→eeγ) | |
| Data 2024 | `/dcache/atlas/mfernand/qt_ntuples/data24/egam4.root` (Z→μμγ) | |
| MC23e | `ntuples/mc23e_egam3_v3/mc23e_egam3_v5.root` | Z→eeγ |
| MC23e | `ntuples/mc23e/mc23e_700771_Zmumug.root` | Z→μμγ |

Channel `llgamma` chains both Z→eeγ and Z→μμγ simultaneously. This is the channel used for all Layer 2 results.

---

## 3. Cluster Geometry

| Parameter | Value |
|-----------|-------|
| Layer | EM Layer 2 |
| Cluster window | 7 eta × 11 phi = **77 cells** |
| Cell size | Δη = 0.025, Δφ = 0.0245 |
| Central cell index | 38 (0-based), = eta=3, phi=5 |
| Branch | `photon.7x11ClusterLr2E` (row-major: phi + 11×eta) |

**Cluster health cut**: require exactly 77 cells present with the central cell being the hottest (highest-energy) cell.

---

## 4. Event Selection

| Cut | Value |
|-----|-------|
| Photon pT | > 10 GeV |
| Photon |η| | < 2.37 (excluding crack 1.37–1.52) |
| Photon ID | None ("no ID") |
| ΔR(γ, ℓ) | > 0.4 |
| m_ℓℓ | 40–83 GeV |
| m_ℓℓγ | 80–100 GeV |
| Pileup reweighting | Applied; MC reweighted to data luminosity |

Scenarios studied:

| Scenario | Conversion | Isolation |
|----------|-----------|-----------|
| unconverted | Unconverted photons only | depends on variant |
| converted | Converted photons only | depends on variant |
| inclusive | Both (hadd of unconverted + converted) | same as components |

---

## 5. Correction Methods

### M1 — Shift Only (flat shift)

For each cell $k$ in the 7×11 cluster:

$$E_k^{\text{MC-RW}} = E_k^{\text{MC}} + \Delta_k \cdot E_{\text{cluster}}^{\text{MC}}$$

where $\Delta_k = \bar{f}_k^{\text{data}} - \bar{f}_k^{\text{MC}}$ (difference of mean normalized cell energy fractions).

### M2 — Shift + Stretch (Francisco's method)

$$E_k^{\text{MC-RW}} = E_{\text{cluster}}^{\text{MC}} \cdot \text{shift}_k + \text{stretch}_k \cdot E_k^{\text{MC}}$$

where:
- $\text{shift}_k = \bar{f}_k^{\text{data}} - \text{stretch}_k \cdot \bar{f}_k^{\text{MC}}$
- $\text{stretch}_k = \sigma_k^{\text{data}} / \sigma_k^{\text{MC}}$ (RMS ratio of normalized energy fraction distributions)

A fallback to M1 ($\text{stretch}_k = 1$) is applied for cells with $\bar{f}_k^{\text{MC}} < 0.002$ (peripheral cells with very low average occupancy).

Final rescaling (both methods): $E_k^{\text{MC-RW}} \mathrel{*}= \sum E_k^{\text{MC}} / \sum E_k^{\text{MC-RW}}$ to guarantee constant cluster energy.

---

## 6. Two-Pass Pipeline (`fill_histograms.C`)

**Pass 1** — Loop data + MC:
- Accumulate per-cell fraction statistics: mean $\bar{f}_k$ and RMS $\sigma_k$ in each (η, pT) bin.
- Fill uncorrected shower shape histograms for data and MC.

**Pass 2** — Loop MC only:
- Compute M1 and M2 correction matrices from Pass 1 statistics.
- Apply corrections cell by cell.
- Recompute $R_\eta$, $R_\phi$, $w_{\eta 2}$ from corrected cells.
- Fill corrected histograms (`mc_M1`, `mc_M2`).

### Histogram Naming Convention

```
h_{var}_{tag}                       # integrated
h_{var}_{tag}_eta{NN}               # per-eta bin
h_{var}_{tag}_eta{NN}_pt{PP}        # per-(eta,pT) bin
cell_profiles/h_frac_{type}_eta{NN} # 7×11 mean cell fraction maps
corrections/h_delta_eta{NN}         # M1 shift map
corrections/h_shift_eta{NN}         # M2 shift map
corrections/h_stretch_eta{NN}       # M2 stretch map
```

Variables: `reta`, `rphi`, `weta2`.
Tags: `data`, `mc`, `mc_M1`, `mc_M2`, `data_stored`, `mc_fudged`, `mc_unfudged`.

---

## 7. Binning Variants

Four variants were run for Layer 2, stored under `output/Layer_2/{variant}/`:

| Variant | Binning | Isolation | Correction dimensions |
|---------|---------|-----------|-----------------------|
| `eta_loose` | 14 eta bins | Loose | [14][1][77] |
| `eta_tight` | 14 eta bins | Tight | [14][1][77] |
| `eta_pt_loose` | 14 eta × 6 pT bins | Loose | [14][6][77] |
| `eta_pt_tight` | 14 eta × 6 pT bins | Tight | [14][6][77] |

**Eta bins** (14 bins, with crack gap at 1.37–1.52):
```
[0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.37, 1.52, 1.6, 1.80, 2.0, 2.2, 2.4]
```

**pT bins** (6 bins, GeV):
```
[10, 15, 20, 25, 30, 40, 1000]
```

---

## 8. Output Structure

For each variant/channel/scenario:

```
output/Layer_2/{variant}/{channel}/{scenario}/
├── histograms.root           # All histograms (Pass 1 + Pass 2)
├── fill.log                  # fill_histograms.C stdout
├── plots/
│   ├── rew_reta.pdf          # Per-eta: Data, MC, M1, M2 (14 pages)
│   ├── rew_rphi.pdf
│   ├── rew_weta2.pdf
│   ├── rew_integrated.pdf    # All-eta integrated (3 pages: reta/rphi/weta2)
│   ├── rew_reta_pt{PP}.pdf   # Per-(eta,pT): 1 PDF per pT bin (eta_pt variants only)
│   ├── rew_rphi_pt{PP}.pdf
│   ├── rew_weta2_pt{PP}.pdf
│   ├── cell_{data,mc,mc_m1,mc_m2}.pdf  # 7×11 heatmaps per eta
│   ├── cell_shift.pdf        # M1 shift map per eta
│   ├── cell_stretch.pdf      # M2 stretch map per eta
│   ├── computed_vs_stored.pdf    # Cell-computed vs branch stored (integrated)
│   ├── computed_vs_stored_eta.pdf
│   ├── fudge_factors.pdf         # Stored branch vs fudged/unfudged MC (integrated)
│   └── fudge_factors_eta.pdf
└── report/
    └── chi2_tables.tex       # χ² per eta bin and shower shape variable
```

Compendium PDFs (generated by `output/Layer_2/make_compendiums.py`):
```
result_compendium_{variant}_{channel}_{scenario}.pdf
```
e.g. `result_compendium_eta_loose_llgamma_unconverted.pdf`

---

## 9. How to Run the Full Layer 2 Pipeline

```bash
# Setup environment
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

# Run everything (all 4 variants × 3 scenarios)
cd scripts/data_mc
bash run_layer2_final.sh

# Or run a subset
VARIANTS="eta_loose" bash run_layer2_final.sh

# Regenerate plots only (skip fill if histograms.root already exists)
OUTBASE="../../output/Layer_2"
for var in eta_loose eta_tight eta_pt_loose eta_pt_tight; do
    for sc in unconverted converted inclusive; do
        root -l -b -q "plot_shower_shapes.C(\"llgamma\", \"${sc}\", \"${OUTBASE}/${var}\", \"${var/_*/}\", \"${var/*_}\")" &
    done
done
wait

# Regenerate compendium PDFs
cd ../../output/Layer_2
python3 make_compendiums.py
```

---

## 10. Shower Shape Variables

| Variable | Formula | Cell window | Branch |
|----------|---------|-------------|--------|
| $R_\eta$ | $E_{3\times 7} / E_{7\times 7}$ | 3×7 numerator, 7×7 denominator | `photon.reta` / `photon.unfudged_reta` |
| $R_\phi$ | $E_{3\times 3} / E_{3\times 7}$ | 3×3 numerator, 3×7 denominator | `photon.rphi` / `photon.unfudged_rphi` |
| $w_{\eta 2}$ | Energy-weighted η width | 3×5 window | `photon.weta2` / `photon.unfudged_weta2` |

These variables are computed both from branches (stored/fudged/unfudged) and directly from the 7×11 cell energy array for cross-validation.

---

## 11. Key Script Summary

| Script | Purpose |
|--------|---------|
| `scripts/data_mc/config.h` | Cluster geometry, eta/pT binning, selection cuts, branch names |
| `scripts/data_mc/fill_histograms.C` | Two-pass pipeline: accumulate cell stats → derive M1/M2 corrections → fill histograms |
| `scripts/data_mc/plot_shower_shapes.C` | Shower shape overlay plots + ratio panels |
| `scripts/data_mc/plot_cell_profiles.C` | 7×11 cell heatmaps + shift/stretch maps |
| `scripts/data_mc/extract_chi2.C` | χ² tables per variant |
| `scripts/data_mc/run_layer2_final.sh` | Full pipeline driver for all 4 variants |
| `output/Layer_2/make_compendiums.py` | LaTeX compendium PDFs from plots |

---

## 12. Known Findings and Decisions

- **M2 outperforms M1** in most η bins: adding the stretch (RMS matching) reduces the MC/Data ratio deviation beyond what shift-only achieves, particularly for $w_{\eta 2}$.
- **Loose isolation variants** were adopted as the primary results. Tight isolation variants are available for systematic studies.
- **η×pT variants** show consistent corrections with η-only variants; the pT-dependence of the correction is mild but visible in the lowest pT bin (10–15 GeV).
- **Converted photons** show larger initial data/MC discrepancies than unconverted, but M2 corrects them similarly well.
- **Crack region** (1.37–1.52 in |η|): excluded from corrections. The crack bin is present in the histogram arrays but has no data/MC content.
- **Peripheral cells** (mean fraction < 0.2%): stretch falls back to 1 (M1 behaviour) to avoid numerical instability.
- **Cell-computed vs. stored shower shapes**: validated that the cell-computed values match the branch-stored values within statistical noise — confirms the cell window definitions are correctly implemented.

---

## 13. What Changes for Layer 1

Layer 1 uses the **strip layer** (EM Layer 1), which has a different cell geometry:

| Parameter | Layer 2 | Layer 1 (expected) |
|-----------|---------|--------------------|
| Cell size Δη | 0.025 | ~0.003125 (8× finer in η) |
| Cell size Δφ | 0.0245 | ~0.098 (4× coarser in φ) |
| Cluster window | 7×11 (η×φ) | different — TBD |
| Branch | `photon.7x11ClusterLr2E` | new branch to be added to NTupleMaker |
| Shower shapes | $R_\eta$, $R_\phi$, $w_{\eta 2}$ | strip-layer shapes: $f_1$, $w_{\eta 1}$, $w_{s,\text{tot}}$, $\Delta E$, $E_\text{ratio}$, $f_\text{side}$ |

### Anticipated Changes to the Pipeline

1. **NTupleMaker** (`cellntuplemaker` GitLab): add new branch for Layer 1 cell energies.
2. **Grid reprocessing**: rerun on MC23e and data24 with updated NTupleMaker to produce ntuples with Layer 1 cells.
3. **`config.h`**: update `kEtaSize`, `kPhiSize`, `kClusterSize`, `kCentralCell`, `kDeltaEta`, `kDeltaPhi` for Layer 1 geometry.
4. **`fill_histograms.C`**: the two-pass logic is geometry-agnostic — only the cluster window constants need updating.
5. **`plot_shower_shapes.C`**: add Layer 1 strip variables ($f_1$, $w_{\eta 1}$, etc.) as new VarDef entries.
6. **Output directory**: use `output/Layer_1/{variant}/` to keep Layer 1 and Layer 2 results separate.
7. **`run_layer2_final.sh`**: clone as `run_layer1_final.sh` with updated output base.
8. **`make_compendiums.py`**: update `OUTBASE` or add a Layer 1 section.

The correction formulas (M1, M2) and the two-pass histogram-filling logic remain unchanged.
