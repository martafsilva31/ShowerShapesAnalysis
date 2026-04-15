# Shower Shapes Analysis — Full Project Context

> **Purpose**: Complete context dump for handoff to other chat sessions or agents.
> **Last updated**: 2026-04-14
> **Author**: Marta Fernandes (mfernand@cern.ch)
> **Workspace**: `/project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis/`

---

## 1. Physics Goal

Derive **data-driven corrections for photon shower shapes** (R_eta, R_phi, w_eta_2) in the ATLAS electromagnetic calorimeter (Layer 2). The corrections fix data-MC disagreements at the **cell energy level** rather than at the shower shape level, so that all observables built from cells are improved simultaneously.

This is more powerful than "fudge factors" (which correct each shower shape variable independently) because a single set of cell-energy corrections automatically improves all downstream variables.

**Reference**: ATL-COM-PHYS-2021-640 (Sections 5.1–5.4).

---

## 2. Environment

| Item | Value |
|------|-------|
| Machine | lxplus (CERN, EL9) |
| ROOT | 6.30.04 via LCG_105a |
| Setup command | `source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh` |
| Working dir | `/project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis/scripts/data_mc/` |
| Output base | `../../output/cell_energy_reweighting_Francisco_method/data24/{channel}/{scenario}/` |
| NTuple source | `/dcache/atlas/mfernand/qt_ntuples/data24/` |
| Git remote | `github.com/martafsilva31/ShowerShapesAnalysis` |

---

## 3. Data & MC Samples

### Input NTuples (produced by NTupleMaker Athena algorithm)

| File | Content | DSID | Events |
|------|---------|------|--------|
| `egam3.root` | Data24 Z→eeγ (EGAM3 derivation) | — | ~4M |
| `egam4.root` | Data24 Z→μμγ (EGAM4 derivation) | — | ~? |
| `mc_eegamma.root` | MC23e Sherpa Z→eeγ | 700770 | ~9.2M |
| `mc_mumugamma.root` | MC23e Sherpa Z→μμγ | 700771 | ~? |

### MC Normalisation

| DSID | Process | σ (pb) | AMI EVNT sumW | totalEvents |
|------|---------|--------|---------------|-------------|
| 700770 | Sh_2214_eegamma | 102.60 | 5.269047e+14 | 185,402,000 |
| 700771 | Sh_2214_mumugamma | 102.59 | 5.421893e+14 | 189,985,000 |

**Luminosity**: 108 fb⁻¹ (data24)

**MC weight formula**: `w = mcwgt × muwgt × (σ × L / sumW_AMI) × SF_l1 × SF_l2 × SF_phIso`

**CRITICAL BUG (now fixed)**: Previously, sumW was read from `h_sumW->GetBinContent(2)` in the NTupleMaker output ROOT files. This histogram only counts DAOD-level events (~36% of the true total after derivation skimming), inflating MC weights by ~2.77× and giving MC/Data ≈ 2.19. The fix uses hardcoded AMI EVNT-level sumW values, giving MC/Data ≈ 0.792. The remaining 20% deficit is due to the incomplete EVNT sample (97.59% coverage, plus filter efficiencies).

---

## 4. Cluster Geometry

- **Layer 2 cluster**: 7 (η) × 11 (φ) = 77 cells
- **Cell size**: Δη = 0.025, Δφ = 0.0245
- **Central cell**: index 38 (ieta=3, iphi=5, 0-based, row-major)
- **Branch**: `photon.7x11ClusterLr2E` — `vector<double>` of 77 values
- **Healthy cluster**: all 77 cells present, central cell has highest energy

---

## 5. Shower Shape Definitions

| Variable | Formula | Cell Window |
|----------|---------|-------------|
| R_eta | E(3×7) / E(7×7) | 3η×7φ inside 7η×7φ |
| R_phi | E(3×3) / E(3×7) | 3η×3φ inside 3η×7φ |
| w_eta_2 | √(Σ E_i × η_i² / Σ E_i − (Σ E_i × η_i / Σ E_i)²) | 3η×5φ window |

**w_eta_2 note**: Our code applies the ATLAS position-dependent polynomial correction (`weta2Correct()` with P0+P1·u+P2·u² coefficients). Francisco's `photoncellbasedrw` does **not** apply this correction — it uses the raw variance with absolute detector eta. This is a known methodological difference (see Section 10).

### NTuple Branches

| Quantity | Branch | Notes |
|----------|--------|-------|
| Stored (fudged for MC) | `photon.reta`, `photon.rphi`, `photon.weta2` | Official ATLAS values |
| Unfudged (MC only) | `photon.unfudged_reta/rphi/weta2` | Before fudge factors |
| Cell energies | `photon.7x11ClusterLr2E` | Layer 2, 77 cells |
| Cluster η for binning | `photon_cluster.eta2` | |
| Cell η for w_eta_2 | `photon_cluster.etamax2` | Used in sub-cell correction |

---

## 6. Event Selection

| Cut | Value | Notes |
|-----|-------|-------|
| Photon pT | > 10 GeV | Francisco's default; was 7 GeV in earlier iterations |
| Photon \|η\| | < 2.37 | Crack excluded [1.37, 1.52) |
| m_ll | [40, 83] GeV | Z window |
| m_llγ | [80, 100] GeV | Radiative Z window |
| ΔR(γ, lepton) | > 0.4 | Both leptons |
| ΔR(l₁, l₂) | No cut | (Francisco's default) |
| Photon ID | No ID requirement | (Francisco's default) |
| Isolation | Loose | (Francisco's default) |
| Cluster health | 77 cells + central hottest | Required |

### Scenarios (selection variants)

| Scenario | Conversion | ID | Isolation |
|----------|-----------|-----|-----------|
| `baseline` | Unconverted only | No ID | Loose |
| `converted` | Converted only | No ID | Loose |
| `all_conv` | Both | No ID | Loose |
| `iso_tight` | Unconverted only | No ID | Tight |

**Channels**: `eegamma` (Z→eeγ), `mumugamma` (Z→μμγ), `llgamma` (both combined)

**Active matrix**: 3 channels × 3 scenarios (baseline, converted, all_conv). `iso_tight` is being dropped.

---

## 7. Correction Methods

### 7.1 Method M1 — Flat Shift (mean-only)

$$E'_k = E_k + \Delta_k \cdot E_{\text{total}}, \quad \Delta_k = \langle f_{\text{data},k} \rangle - \langle f_{\text{MC},k} \rangle$$

where $f_k = E_k / E_{\text{total}}$ is the normalised cell energy fraction.

- Corrects only the mean of each cell's energy distribution
- Preserves total energy exactly (Σ Δ_k = 0 by construction)
- No rescaling needed

### 7.2 Method M2 — Shift + Stretch (mean + variance)

$$E'_k = E_{\text{total}} \cdot \text{shift}_k + \text{stretch}_k \cdot E_k$$

where:
- $\text{stretch}_k = \sigma_{\text{data},k} / \sigma_{\text{MC},k}$
- $\text{shift}_k = \mu_{\text{data},k} - \text{stretch}_k \cdot \mu_{\text{MC},k}$

After application, energy conservation rescaling:
$$E'_k \mathrel{*}= \frac{E_{\text{total}}}{\sum_k E'_k}$$

**This is mathematically identical to Francisco's method** (confirmed by code comparison). The equivalence in normalised form:

$$f_k^{\text{rw}} = \frac{\sigma_{\text{data},k}}{\sigma_{\text{MC},k}} (f_k^{\text{MC}} - \mu_{\text{MC},k}) + \mu_{\text{data},k}$$

### Correction Derivation (Pass 1)

Statistics accumulated per eta bin, per cell:
- **Data**: $\sum f_k$, $\sum f_k^2$, count
- **MC**: $\sum w \cdot f_k$, $\sum w \cdot f_k^2$, $\sum w$ (MC-weighted)
- Mean: $\mu = \overline{f_k}$, Variance: $\sigma^2 = \overline{f_k^2} - \overline{f_k}^2$

Corrections stored per eta bin per cell: `delta[eta][cell]` (M1), `shift[eta][cell]`, `stretch[eta][cell]` (M2).

### Binning for Corrections

**Our code**: 14 |η| bins: {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.37, 1.52, 1.6, 1.80, 2.0, 2.2, 2.4}
- Bin 8 [1.37, 1.52) is the crack region — always 0 events, corrections skipped.

**Francisco's code** (`photoncellbasedrw`): 6 pT bins × 5 hit-position bins = 30 bins per conversion type. NO eta binning.
- pT: [10, 15, 20, 25, 30, 40, 1000] GeV
- Hit position: 1 central region (λ < threshold) + 4 quadrants (sign of Δη × sign of Δφ)

These are **completely orthogonal** decompositions addressing different physics:
- η binning captures detector geometry variations
- pT binning captures energy-scale effects
- Hit-position binning captures impact-point effects on cell response

---

## 8. Code Architecture

### Files in `scripts/data_mc/`

| File | Lines | Purpose |
|------|-------|---------|
| `config.h` | 479 | Shared config: geometry, selection, branches, shower shape formulas, MC weights |
| `fill_histograms.C` | 614 | Two-pass pipeline: Pass 1 (accumulate stats) → compute corrections → Pass 2 (apply corrections) |
| `plot_shower_shapes.C` | 375 | Data vs MC shower shape comparison plots with ratio panel. SET A: branch (fudged). SET B: cell-computed (M1/M2). Supports nRebin parameter. |
| `plot_kinematics.C` | 285 | pT, η distributions + zero-padding diagnostic plots |
| `plot_cell_profiles.C` | 318 | 7×11 cell fraction heatmaps: data mean, MC mean, M1-corrected, M2-corrected, shift, stretch |
| `extract_chi2.C` | 367 | χ² comparison of data vs MC/M1/M2 distributions |
| `make_compendium.py` | 398 | LaTeX compendium PDF per scenario (collects all plots) |
| `run.sh` | 110 | Driver: fill + plot pipeline, supports `--batch` and `--plot-only` |

### Pipeline Flow

```
                     fill_histograms.C
                    ┌───────────────────────────────────┐
                    │  Pass 1: Data + MC                │
                    │  → accumulate cell fraction stats  │
                    │  → fill uncorrected SS histos      │
                    │  → compute M1 + M2 corrections     │
                    │                                   │
                    │  Pass 2: MC only                   │
                    │  → apply M1, M2 corrections        │
                    │  → fill corrected SS histos        │
                    │  → build cell profile heatmaps     │
                    └────────────┬──────────────────────┘
                                 │
                                 ▼
                       histograms.root
                                 │
              ┌──────────────────┼──────────────────┐
              ▼                  ▼                  ▼
   plot_shower_shapes.C   plot_kinematics.C   plot_cell_profiles.C
         │                       │                  │
         ▼                       ▼                  ▼
    rew_*.pdf              kinematics_*.pdf    cell_*.pdf
    fudge_factors.pdf      zero_padding.pdf
    computed_vs_stored.pdf
```

### Driver Usage

```bash
# Single scenario
./run.sh eegamma baseline
./run.sh eegamma baseline --plot-only

# All 3 channels × 3 scenarios
./run.sh --batch
./run.sh --batch --plot-only

# Direct ROOT invocation
root -l -b -q 'fill_histograms.C("eegamma", "baseline")'
root -l -b -q 'plot_shower_shapes.C("eegamma", "baseline", "../../output/cell_energy_reweighting_Francisco_method/data24")'
root -l -b -q 'plot_shower_shapes.C("eegamma", "baseline", "../../output/cell_energy_reweighting_Francisco_method/data24", 2)'  # nRebin=2 → 50 bins
root -l -b -q 'plot_shower_shapes.C("eegamma", "baseline", "../../output/cell_energy_reweighting_Francisco_method/data24", 4)'  # nRebin=4 → 25 bins
```

### Histogram Naming Convention

In `histograms.root`:
- **Integrated**: `h_{var}_{tag}` — e.g. `h_reta_data_computed`, `h_reta_mc_M2`
- **Per-eta bin**: `h_{var}_{tag}_eta{NN}` — e.g. `h_reta_data_eta00`, `h_reta_mc_M1_eta07`
- **Cell profiles** (in `cell_profiles/` subdir): `h_frac_{type}_eta{NN}`

Variable names: `reta`, `rphi`, `weta2`
Tags: `data`, `data_stored`, `mc` (cell-computed), `mc_fudged`, `mc_unfudged`, `mc_M1`, `mc_M2`

---

## 9. Current Configuration (as of 2026-04-14)

| Parameter | Value | Notes |
|-----------|-------|-------|
| `kNBins` | 100 | Changed from 50 on 2026-04-14 |
| `kExcludeZeroPadFromCorr` | `true` | Added 2026-04-14. Excludes events with zero-padded cells from correction derivation, but still includes them in shower shape histograms. |
| `kClusterSize` | 77 | 7×11 |
| `kCentralCell` | 38 | (ieta=3, iphi=5) |
| `kNEtaBins` | 14 | Crack at bin 8 |
| nRebin options | 1, 2, 4 | Produces 100-bin, 50-bin, 25-bin plots |

### Zero-Padding

- Data: 0.01% events have ≥1 zero cell (4 / 78,131 for eegamma/baseline)
- MC: 1.34% events have ≥1 zero cell (12,818 / 955,581 for eegamma/baseline)
- With `kExcludeZeroPadFromCorr = true`, these events are excluded from Pass 1 correction derivation but still contribute to shower shape histograms.

### Output Status

All 12 scenarios (3 channels × 4) were filled on 2026-04-13 with 50-bin config.
**Only eegamma/baseline** has been refilled on 2026-04-14 with the new 100-bin + zero-pad exclusion config.
The remaining 8 active scenarios (3ch × 3sc minus eegamma/baseline) still use old 50-bin config.

---

## 10. Comparison with Francisco's `photoncellbasedrw`

Full comparison document: [`docs/francisco_comparison.md`](francisco_comparison.md)

### Summary of Matches ✅

1. **Core M2 formula** = Francisco's shift+stretch (mathematically identical)
2. **Cell grid**: 7×11, central cell 38, health check — identical
3. **Mass windows**: m_ll ∈ [40, 83] GeV, m_llγ ∈ [80, 100] GeV
4. **ΔR cut**: > 0.4
5. **Energy conservation rescaling**: Both rescale after M2
6. **Cell fraction normalisation**: E_k / E_total for all 77 cells
7. **MC weighting** in profile accumulation

### Key Differences ⚠️

| Aspect | Our code | Francisco's `photoncellbasedrw` | Impact |
|--------|----------|----------------------------------|--------|
| Histogram bins | 100 (was 50) | 100 | Now matched |
| Correction binning | 14 η bins | 6 pT × 5 hitPos = 30 bins | **HIGH** — completely orthogonal decompositions |
| w_eta_2 calculation | Relative η + ATLAS polynomial correction | Absolute detector η, no correction | **MEDIUM** — different w_eta_2 values |
| Zero-pad handling | Events with zero cells excluded from corrections | Not explicitly addressed | **LOW** |
| pT threshold | 10 GeV | 10 GeV | Matched |
| Photon ID | No ID (baseline) | No ID | Matched |
| Isolation | Loose (baseline) | Loose | Matched |

### Reference Code Locations

| Code | Path | Description |
|------|------|-------------|
| `photoncellbasedrw` | `/project/atlas/users/mfernand/QT/ShowerShapes/photoncellbasedrw/` | Francisco's clean Python+C++ framework |
| `showershapereweighting` | `/project/atlas/users/mfernand/QT/ShowerShapes/showershapereweighting/` | Francisco's older C++ code (CalCoef.C, NTUP.C), uses 14 eta bins |
| `NTupleMaker` | `/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace/source/NTupleMaker/` | Athena C++ algorithm producing cell-level ntuples |

---

## 11. Known Issues & Diagnostics

### Active Issues

1. **Reweighting not performing as expected**: M1 and M2 corrections improve Data-MC agreement somewhat, but residual differences remain, especially in w_eta_2. Root causes identified:
   - **No pT binning** (HIGH impact): Our corrections are derived in η bins only. Francisco's code uses pT × hit-position bins. The absence of pT dependence means corrections are averaged over a wide energy range.
   - **w_eta_2 polynomial correction** (MEDIUM impact): We apply the ATLAS position correction; Francisco does not. This creates a systematic difference in w_eta_2 that may propagate to correction quality.
   - **Zero-padding asymmetry** (LOW impact, now mitigated): MC has 100× more zero-padded cells than data. Now excluded from correction derivation.

2. **MC/Data normalisation = 0.792**: ~20% deficit after AMI sumW fix. Caused by incomplete EVNT sample (97.59% of total events generated), filter efficiencies, and potentially missing MC campaigns. Not a bug — the shape comparisons are unaffected (histograms normalised to area for shape comparisons).

3. **Data22 NTuples have 0 events**: Runs 428648–429027 (periods A–E) are not in the GRL which starts at run 431810 (period F). Would need to resubmit for data22 periods F/H/J. Currently using data24 only.

### Resolved Issues

1. **sumW from h_sumW inflating MC weights by 2.77×** (FIXED 2026-04-13): Replaced GetBinContent(2) from DAOD-level h_sumW with hardcoded AMI EVNT-level sumW values.

---

## 12. Planned Work (Overnight Diagnostic Plan)

Organised into phases with dependencies. The plan was agreed on 2026-04-13.

### Phase A: Code Fixes (in progress)

| Step | Description | Status |
|------|-------------|--------|
| A1 | Increase kNBins from 50 → 100 | ✅ Done |
| A2 | Add kExcludeZeroPadFromCorr flag | ✅ Done |
| A3 | Remove w_eta_2 ATLAS polynomial correction (match Francisco) | ❌ Pending |

### Phase B: New Diagnostic Plot Scripts

| Script | Purpose | Status |
|--------|---------|--------|
| Fudge factors per-eta | Plot stored (fudged) vs unfudged vs data per eta bin | ❌ Pending |
| Computed vs stored per-eta | Verify cell-computed SS match stored branches per eta | ❌ Pending |
| Chi-squared table | Quantitative comparison: χ²(Data, MC/M1/M2) per eta | ❌ Pending |
| Shift/stretch summary | Bar chart of M2 shift and stretch values across eta/cells | ❌ Pending |
| Normalised kinematics | Shape-normalise pT/η distributions for Data vs MC comparison | ❌ Pending |

### Phase C: Refill All Scenarios

Refill all 9 active scenarios (3 channels × 3 scenarios) with 100-bin + zero-pad exclusion config.
Only eegamma/baseline done so far.

### Phase D: Plot All Scenarios

Run all plot scripts on all 9 scenarios.

### Phase E: Compendia & Report

Generate LaTeX compendia (one per scenario) and update the data-MC comparison report.

---

## 13. Typical Event Counts (eegamma/baseline, 2026-04-14 fill)

| Quantity | Value |
|----------|-------|
| Data total entries | 4,018,791 |
| Data passing selection | 78,131 |
| MC total entries | 9,249,301 |
| MC passing selection | 955,581 |
| Data zero-padded | 4 (0.01%) |
| MC zero-padded | 12,818 (1.34%) |

### Events per η bin (data / MC weighted)

| Bin | Range | Data | MC (weighted) |
|-----|-------|------|---------------|
| 0 | [0.00, 0.20) | 9,813 | 7,827 |
| 1 | [0.20, 0.40) | 11,648 | 9,365 |
| 2 | [0.40, 0.60) | 11,245 | 9,114 |
| 3 | [0.60, 0.80) | 10,847 | 8,727 |
| 4 | [0.80, 1.00) | 9,089 | 7,146 |
| 5 | [1.00, 1.20) | 7,450 | 5,805 |
| 6 | [1.20, 1.30) | 3,222 | 2,577 |
| 7 | [1.30, 1.37) | 704 | 648 |
| 8 | [1.37, 1.52) | 0 | 0 |
| 9 | [1.52, 1.60) | 1,911 | 1,153 |
| 10 | [1.60, 1.80) | 3,873 | 2,705 |
| 11 | [1.80, 2.00) | 3,333 | 2,472 |
| 12 | [2.00, 2.20) | 3,013 | 2,251 |
| 13 | [2.20, 2.40) | 1,979 | 1,344 |

---

## 14. Output Files per Scenario

After a full `run.sh` invocation, each `{channel}/{scenario}/` directory contains:

```
histograms.root           # All histograms (integrated + per-eta + cell profiles)
plots/
  rew_reta.pdf            # SET B: R_eta (Data cell vs MC cell vs M1 vs M2), 100 bins
  rew_reta_50bins.pdf     # Same, rebinned to 50 bins
  rew_reta_25bins.pdf     # Same, rebinned to 25 bins
  rew_rphi.pdf            # SET B: R_phi, 100 bins
  rew_rphi_50bins.pdf     # 50 bins
  rew_rphi_25bins.pdf     # 25 bins
  rew_weta2.pdf           # SET B: w_eta_2, 100 bins
  rew_weta2_50bins.pdf    # 50 bins
  rew_weta2_25bins.pdf    # 25 bins
  fudge_factors.pdf       # Fudge factor comparison, 100 bins
  fudge_factors_50bins.pdf
  fudge_factors_25bins.pdf
  computed_vs_stored.pdf  # Cell-computed vs branch-stored SS, 100 bins
  computed_vs_stored_50bins.pdf
  computed_vs_stored_25bins.pdf
  cell_data.pdf           # 7×11 heatmap: data mean cell fraction (14 pages, one per η bin)
  cell_mc.pdf             # 7×11 heatmap: MC mean cell fraction
  cell_mc_m1.pdf          # 7×11 heatmap: MC after M1 correction
  cell_mc_m2.pdf          # 7×11 heatmap: MC after M2 correction
  cell_shift.pdf          # 7×11 heatmap: M2 shift values
  cell_stretch.pdf        # 7×11 heatmap: M2 stretch values
  kinematics_pt.pdf       # Photon pT distribution (Data vs MC)
  kinematics_eta.pdf      # Photon |η| distribution (Data vs MC)
  zero_padding.pdf        # Number of zero cells per event (Data vs MC)
```

---

## 15. Key Formulas Quick Reference

### MC event weight
```
w = mcwgt × muwgt × (σ_DSID × L / sumW_DSID) × SF_l1 × SF_l2 × SF_phIso
```

### Cell fraction
```
f_k = E_k / Σ_k E_k      (k = 0..76, normalised to 77-cell total)
```

### M1 correction (flat shift)
```
E'_k = E_k + Δ_k × E_total
Δ_k = <f_data,k> - <f_MC,k>
```

### M2 correction (shift + stretch)
```
stretch_k = σ_data,k / σ_MC,k
shift_k = μ_data,k - stretch_k × μ_MC,k
E'_k = E_total × shift_k + stretch_k × E_k
E'_k *= E_total / Σ E'_k    (energy conservation)
```

### Shower shapes from cells
```
R_eta = E(3×7) / E(7×7)
R_phi = E(3×3) / E(3×7)
w_eta_2 = sqrt(Σ E_i × η_i²/Σ E_i - (Σ E_i × η_i/Σ E_i)²)   [3×5 window]
```

### Healthy cluster check
```
size == 77  AND  E_total > 0  AND  E_38 ≥ E_k for all k ≠ 38
```

---

## 16. Git Status & Conventions

- **Branch**: `master`
- **Commit style**: Conventional commits (`feat:`, `fix:`, `docs:`, `chore:`)
- **Never commit**: `.root`, `ntuples/`, `output/`, `*.log`, `*.pdf`, `*.png`, `*.pyc`
- **Recent uncommitted changes**:
  - `fill_histograms.C`: kNBins 50→100, kExcludeZeroPadFromCorr flag
  - `run.sh`: Added nRebin=4 variant

---

## 17. Sibling Repositories

| Repo | Path | Remote | Purpose |
|------|------|--------|---------|
| NTupleMaker | `../NTupleMaker_workspace/source/NTupleMaker/` | `gitlab.cern.ch/femarta/cellntuplemaker` | Athena C++ algorithm producing cell-level ntuples from xAOD |
| showershapereweighting | `../showershapereweighting/` | `gitlab.cern.ch/mobelfki/showershapereweighting` | Supervisor's (Francisco's) Run 2 C++ correction code |
| photoncellbasedrw | `../photoncellbasedrw/` | (Francisco's) | Clean Python+C++ framework — the primary reference implementation |

---

## 18. Instructions for Continuing This Work

1. **Source the environment** before any ROOT command:
   ```bash
   source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
   ```

2. **Work from the scripts directory**:
   ```bash
   cd /project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis/scripts/data_mc/
   ```

3. **Fill time**: ~14 minutes per scenario for eegamma (4M data + 9.2M MC entries).

4. **Review the analysis conventions** in `.github/` before modifying code.

5. **Immediate next steps** (as of 2026-04-14):
   - Review eegamma/baseline plots at 100/50/25 bins to decide preferred binning
   - Implement Phase A3 (remove w_eta_2 polynomial correction)
   - Refill remaining 8 scenarios with new config
   - Implement Phase B diagnostic plots
   - Generate compendia and update LaTeX report
