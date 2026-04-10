# Comparison: `photoncellbasedrw` (Francisco) vs `ShowerShapesAnalysis` (Marta)

## Executive Summary

Francisco's repository implements **one method** — a shift+stretch correction on normalised cell energies, using pT × hit-position binning with separate unconverted/converted treatments.
**This method is mathematically identical to your M3 (shift+stretch).**
Your M1 (flat shift) and M2 (ΔR-matched TProfile) have no counterpart in Francisco's code.

The key differences are the **binning axes** (Francisco: pT × hit-position × conversion; you: η only, unconverted only) and the **shower shape computation** (Francisco: grid-index windows; you: position-based windows with an ATLAS w_eta2 sub-cell correction).

---

## 1. Repository Structure

| Aspect | `photoncellbasedrw` (Francisco) | `ShowerShapesAnalysis` (Marta) |
|--------|----------------------------------|-------------------------------|
| Language | Python + C++ via `gInterpreter` | C++ ROOT macros |
| Entry point | `PhotonRW.py -c config.yaml` | Individual `.C` scripts / `run_pipeline.sh` |
| Config | Single `config.yaml` (years, binning, method, plots) | `config.h` (constants, cell layout, η binning) |
| Steps | profileCalc → rwCalc → profileRwCalc → plotCalc | `derive_corrections.C` → `validate_data_mc.C` → plotting scripts |
| Old code | `src/old/` (deprecated `normalizedAvg` / `absoluteAvg` methods) | N/A |

---

## 2. Method Comparison

### 2.1 Francisco's Method 1 ≡ Your M3

**Francisco's README formula:**

$$E_{i,j}^{RW} = \frac{\text{RMS}^{\text{data}}_{i}}{\text{RMS}^{\text{MC}}_{i}}\; E_{i,j}^{MC} \;+\; E^{MC}\!\left(\bar{e}_i^{\text{data}} - \frac{\text{RMS}^{\text{data}}_{i}}{\text{RMS}^{\text{MC}}_{i}}\;\bar{e}_i^{\text{MC}}\right)$$

**In Francisco's code** (`RWCalculator.py`):
```python
# correctionFunction():
stretch = rms_e_data / rms_e_mc        # ROOT GetRMS() = std dev (σ)
shift   = e_data - stretch * e_mc      # e_data, e_mc = <f_k> from TProfile

# rw_func():
E_rw = sumE * shift + stretch * cell_E
```

Converting to normalised fractions:

$$f_k^{\text{rw}} = \frac{\sigma_{\text{data},k}}{\sigma_{\text{MC},k}} \left(f_k^{\text{MC}} - \mu_{\text{MC},k}\right) + \mu_{\text{data},k}$$

**Your M3** (`derive_corrections.C` / `validate_data_mc.C`):

$$f'_k = \frac{\sigma_{\text{data},k}}{\sigma_{\text{MC},k}} \left(f_{\text{MC},k} - \mu_{\text{MC},k}\right) + \mu_{\text{data},k}, \qquad E'_k = f'_k \cdot E_{\text{total}}$$

**These are identical.**

> **Note on ROOT `GetRMS()`:** Despite the misleading name, ROOT's `TH1::GetRMS()` returns the sample **standard deviation** (σ), not the root-mean-square. Francisco's `rms_e_data` is actually σ_data, matching your σ computation via √(⟨f²⟩ − ⟨f⟩²).

### 2.2 Your M1 — No counterpart in Francisco's code

Your M1 uses only the **mean shift** (no stretch):

$$E_{\text{new},k} = E_{\text{old},k} + \Delta_k \cdot E_{\text{total}}, \qquad \Delta_k = \langle f_{\text{data},k}\rangle - \langle f_{\text{MC},k}\rangle$$

This is like setting `stretch = 1` in Francisco's formula (no variance matching). M1 preserves total energy exactly (∑Δ_k = 0 by construction), so no rescaling is needed.

Francisco's deprecated `src/old/corrections/corrector.py` had a `normalizedAvg` mode that also did **not** rescale, but the formula details are different (it stored a single correction coefficient per cell rather than separate shift/stretch).

### 2.3 Your M2 — No counterpart in Francisco's code

Your M2 is an entirely different strategy:
- ΔR-match individual data/MC photons (ΔR < 0.1)
- Build a TProfile of α = f_data,k − f_MC,k as a function of f_MC,k per cell per η bin
- At application time, interpolate: α(f_MC,k) → apply E'_k = E_k + α · E_total
- Fall back to M1 if TProfile has < 10 entries per bin

This event-level matching approach is not present anywhere in Francisco's code.

---

## 3. Binning

This is the **largest operational difference** between the two codebases.

| Axis | Francisco | Marta |
|------|-----------|-------|
| **η (pseudorapidity)** | ❌ Not binned | ✅ 14 bins (0.00–2.40, crack at bin 8) |
| **pT (transverse momentum)** | ✅ 6 bins: 10, 15, 20, 25, 30, 40, 1000 GeV | ❌ Not binned |
| **Hit position** | ✅ 5 regions (central + 4 quadrants by sign of Δη/Δφ) | ❌ Not binned |
| **Conversion** | ✅ Separate unconverted (u) / converted (c) | ✅ Unconverted only |
| **Application binning** | Can differ from derivation (pT-inclusive, hitPos-inclusive) | Same as derivation |

**Key implication:** Francisco derives corrections in fine pT bins + hit-position regions, then can apply them inclusively. You derive corrections in fine η bins with no pT dependence. These approaches address orthogonal sources of variation:
- pT binning captures energy-scale dependent effects
- η binning captures detector-geometry dependent effects
- Hit-position binning captures photon-impact-point effects on the cluster shape

---

## 4. Cell Grid & Indexing

Both use the same **7 η × 11 φ = 77 cell** EM Layer 2 cluster. Central cell is index **38** in both.

| Aspect | Francisco | Marta |
|--------|-----------|-------|
| Index convention | `cell = eta * 11 + phi` (0-based), stored as 1-based in histograms | `k = phi + 11 * eta` (0-based) |
| Central cell | 38 → (η=3, φ=5) | 38 → (η=3, φ=5) ✅ Same |
| Health check | `n_L2 == 77` AND max-energy cell == 38 | Same |

**⚠️ Index ordering note:** Francisco iterates `for ieta ... for iphi` and uses `cell = (ieta-1)*11 + iphi - 1` (1-based loop indices). You use `k = phi + kPhiSize * eta` (0-based). Both resolve to the same cell for the same (η,φ) position. No conflict.

---

## 5. Shower Shape Computation

| Variable | Francisco (grid-index windows) | Marta (position-based windows) |
|----------|-------------------------------|--------------------------------|
| **E(3×3)** | η ∈ [2,4], φ ∈ [4,6] (hardcoded) | Geometric: \|Δη\| < 1.5 Δη_cell, \|Δφ\| < 1.5 Δφ_cell relative to hottest cell |
| **E(3×7)** | η ∈ [2,4], φ ∈ [2,8] | Geometric: \|Δη\| < 1.5 Δη_cell, \|Δφ\| < 3.5 Δφ_cell |
| **E(7×7)** | η ∈ [0,6], φ ∈ [2,8] | Geometric: \|Δη\| < 3.5 Δη_cell, \|Δφ\| < 3.5 Δφ_cell |
| **R_η** | E(3×7) / E(7×7) | Same formula |
| **R_φ** | E(3×3) / E(3×7) | Same formula |
| **w_eta2** | 3×5 window (φ [3,7], η [2,4]); no sub-cell correction | 3×5 window (position-based); **ATLAS sub-cell polynomial correction applied** |

**Key difference:** Your w_eta2 includes the standard ATLAS position correction (`P0 + P1·u + P2·u²` where u is the sub-cell hit position), which Francisco's code does not apply. This means your w_eta2 values will differ slightly from Francisco's — not a bug, but a methodological choice.

For E(3×3), E(3×7), E(7×7) the results should be **numerically identical** when the hottest cell is at index 38, since the geometric windows and the hardcoded index windows select the same cells in a regular grid centred on (3,5).

---

## 6. Event Selection

| Cut | Francisco | Marta |
|-----|-----------|-------|
| Photon pT | > 10 GeV | > 7 GeV |
| Photon \|η\| | < 2.37, crack excluded | < 2.37, crack excluded ✅ |
| m_ll | [40, 83] GeV | [40, 83] GeV ✅ |
| m_llγ | [80, 100] GeV | [80, 100] GeV ✅ |
| ΔR(γ, lepton) | > 0.4 | > 0.4 ✅ |
| ΔR(l₁, l₂) | Not applied | > 0.2 |
| Identification | Configurable (default: noID) | Tight required |
| Isolation | Configurable (default: loose) | None |
| Conversion | Separate u/c treatment | Unconverted only |
| MC truth match | Not mentioned | Required |
| Healthy cluster | 77 cells + central hottest | 77 cells + central hottest ✅ |
| Channel | Both (Z→eeγ + Z→μμγ, via ntuples) | Both (separate EGAM3/EGAM4) |

Notable differences:
- **pT threshold**: Francisco is 3 GeV higher (10 vs 7 GeV), which reduces the low-pT noise in corrections
- **ID/Isolation**: Francisco runs with no ID + loose iso; you require tight ID and no isolation cut — these select different photon populations
- **Converted photons**: Francisco processes them separately; you exclude them entirely

---

## 7. Energy Conservation (Rescaling)

| | Francisco | Your M1 | Your M2 | Your M3 |
|--|-----------|---------|---------|---------|
| Rescaling? | **Yes**: `r = sumE / sumE_rw`, then `E_k *= r` | **No** (exact by construction) | **Yes** | **Yes** |
| When skipped | Never (always rescales) | Always (not needed) | — | — |

Francisco **always** rescales after applying the shift+stretch correction, even though the formula should approximately conserve energy. Your M3 also rescales. Match is exact here.

For your M1, the mean-shift formulation guarantees ∑Δ_k = 0, so the total energy is preserved automatically.

---

## 8. Summary of Matches and Mismatches

### ✅ Things that match
1. **Core method**: Francisco's Method 1 = your M3 (shift+stretch), identical formula
2. **Cell grid**: Same 7×11 layout, same central cell (38), same health check
3. **Mass windows**: m_ll ∈ [40, 83] GeV, m_llγ ∈ [80, 100] GeV
4. **ΔR(γ,l) cut**: > 0.4
5. **Rescaling strategy**: Both rescale E_total after M3/shift+stretch application
6. **E(3×3), E(3×7), E(7×7) computation**: Equivalent when cluster is centred on cell 38

### ⚠️ Significant differences
1. **Number of methods**: Francisco has 1 method; you have 3 (M1, M2, M3)
2. **Binning axes**: pT × hit-position (Francisco) vs η (you) — completely orthogonal
3. **Conversion treatment**: Both u/c (Francisco) vs unconverted-only (you)
4. **w_eta2 sub-cell correction**: You apply it; Francisco does not
5. **MC weights in correction derivation**: Francisco uses event weights (pileup, gen, xs, trigger, lumi) during profile accumulation; your M1/M3 use weight = 1.0 for population means
6. **pT threshold**: 10 GeV (Francisco) vs 7 GeV (you)
7. **Photon ID**: noID (Francisco default) vs tight (you)
8. **Isolation**: loose (Francisco default) vs none (you)

### 🔍 Minor differences
1. **Language**: Python vs C++
2. **Input format**: Two separate trees matched by event number (Francisco) vs single tree with cell branches (you)
3. **M3 sigma clamping**: Your M3 clamps σ_ratio to [0.5, 2.0]; Francisco does not clamp
4. **Old deprecated methods**: Francisco's `src/old/` had normalizedAvg and absoluteAvg approaches (no longer used)

---

## 9. Recommendations

1. **Adopt η binning into Francisco's framework** (or vice versa): The combination of pT + η + hit-position binning could give better corrections than either alone.
2. **Align pT threshold and ID/isolation**: Decide on a common selection to ensure the correction population is consistent.
3. **Investigate w_eta2 correction**: Determine whether the ATLAS sub-cell correction should be applied (it's the official recipe; Francisco may have omitted it intentionally or not).
4. **Consider MC weight treatment**: Francisco weights MC events during profile accumulation; you use unit weights. The weighted approach is more correct for drawing physics conclusions, but unit weights give the actual MC distribution "as simulated."
5. **Converted photons**: If the analysis scope includes converted photons, your M1/M2/M3 would need extension.
