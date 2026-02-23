# MC23e Cell-Energy Reweighting Closure Test

Closure test pipeline for the photon shower shape cell-energy reweighting method,
applied to MC23e `Z→eeγ` ntuples. Validates the method described in
[ATL-COM-PHYS-2021-640](https://cds.cern.ch/record/2774785), Section 5.

The pipeline creates **pseudo-data** from MC by applying known per-cell distortions,
derives reweighting corrections, applies them, and checks whether the corrected
MC recovers the pseudo-data distributions. Two correction methods are compared.

## Method

### 1. Pseudo-data generation (`create_pseudodata.C`)

Each cell energy $e_k$ (for cell $k$ in the 7×11 Layer-2 cluster) is distorted:

$$
e'_k = e_k \cdot \left(1 + \delta \cdot b_k + \mathcal{G}(0,\, \sigma)\right)
$$

where:
- $\delta$ is the **systematic distortion level** (e.g. 0.01, 0.03, 0.05)
- $b_k$ is a **fixed bias pattern** based on Chebyshev distance from the central cell:
  | Distance | Bias $b_k$ | Region |
  |----------|-----------|--------|
  | 0        | +1.0      | Central cell |
  | 1        | +0.4      | 1st ring |
  | 2        | −0.4      | 2nd ring |
  | ≥3       | −1.0      | Outer cells |
- $\mathcal{G}(0, \sigma)$ is **per-event per-cell Gaussian noise** with $\sigma = \delta/2$

After scaling, total cluster energy is conserved:

$$
e'_k \leftarrow e'_k \cdot \frac{E_\text{tot}}{\sum_j e'_j}
$$

This creates narrower showers (enhanced centre, suppressed edges), mimicking
typical data–MC differences observed in ATLAS calorimetry.

### 2. Correction derivation (`derive_corrections.C`)

For each event, cell energy fractions are computed:

$$
f_k = \frac{e_k}{E_\text{tot}} = \frac{e_k}{\sum_j e_j}
$$

**Method 1 — Flat shift (§5.1):**

A single additive correction per cell per $|\eta|$ bin:

$$
\Delta_k = \langle f_k^\text{pseudo} \rangle - \langle f_k^\text{MC} \rangle
$$

**Method 2 — Energy-dependent TProfile correction (§5.2):**

A `TProfile` records the mean difference as a function of the MC fraction:

$$
\alpha_k(f) = \langle f_k^\text{pseudo} - f_k^\text{MC} \rangle \Big|_{f_k^\text{MC} = f}
$$

This captures the energy-dependent (linear) relationship between the correction
and the cell fraction, following the approach in the supervisor's `tree.C`.
The TProfile uses the same binning as the supervisor's `var.h`:

```
EBins = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12.5, 15}
```

### 3. Correction application (`apply_corrections.C`)

**Method 1 (Shift Only):**

$$
f_k^\text{corr} = f_k^\text{MC} + \Delta_k
$$

**Method 2 (Shift + Stretch):**

$$
f_k^\text{corr} = f_k^\text{MC} + \alpha_k\!\left(f_k^\text{MC}\right)
$$

where $\alpha_k$ is obtained via `TProfile::Interpolate()`, exactly as in the
supervisor's `tree.C` (line ~197).

In both cases, cluster energy is conserved by rescaling:

$$
e_k^\text{corr} = \frac{f_k^\text{corr}}{\sum_j f_j^\text{corr}} \cdot E_\text{tot}
$$

### 4. Validation (`validate_closure.C`)

Shower shape variables are recomputed from the corrected cell energies:

| Variable | Formula | Window |
|----------|---------|--------|
| $R_\eta$ | $E(3\times7) / E(7\times7)$ | 3η × 7φ / 7η × 7φ |
| $R_\phi$ | $E(3\times3) / E(3\times7)$ | 3η × 3φ / 3η × 7φ |
| $w_{\eta\,2}$ | $\sqrt{\langle\eta^2\rangle_E - \langle\eta\rangle_E^2}$ | 3η × 5φ |

Normalised distributions are compared using $\chi^2$/ndf between the
reweighted MC and pseudo-data, binned in 14 $|\eta|$ bins matching the
supervisor's `var.h`.

### 5. Plotting (`plot_closure.C`)

Multi-page PDFs with:
- **Upper panel**: Overlay of pseudo-data (markers), original MC (red),
  Method 1 (green dashed), Method 2 (blue solid)
- **Lower panel**: Ratio to pseudo-data

One page per $|\eta|$ bin, one PDF per variable.

## File Structure

| File | Purpose |
|------|---------|
| `config.h` | Shared constants, geometry, binning, shower shape functions |
| `create_pseudodata.C` | Generate pseudo-data with systematic + noise distortion |
| `derive_corrections.C` | Derive M1 flat shifts and M2 TProfile corrections |
| `apply_corrections.C` | Apply corrections to MC cell energies |
| `validate_closure.C` | Fill shower shape histograms, compute χ² metrics |
| `plot_closure.C` | Generate comparison plots with ratio panels |
| `run_closure_test.sh` | Single-level driver script (6 pipeline steps) |
| `run_closure_suite.sh` | Multi-level wrapper (runs 1%, 3%, 5%) |

## Running

### Prerequisites

```bash
# On lxplus (CERN)
setupATLAS
lsetup "root 6.32.08-x86_64-el9-gcc13-opt"
```

### Single distortion level

```bash
cd scripts/closure_test/
./run_closure_test.sh 500000 0.05    # 500k events, 5% distortion
```

### Full suite (1%, 3%, 5%)

```bash
./run_closure_suite.sh 500000
```

### Output

```
output/closure_test_1pct/    # ROOT files for 1% level
output/closure_test_3pct/    # ROOT files for 3% level
output/closure_test_5pct/    # ROOT files for 5% level
plots/closure_test_1pct/     # PDFs for 1% level
plots/closure_test_3pct/     # PDFs for 3% level
plots/closure_test_5pct/     # PDFs for 5% level
```

## Design Notes

### Matching the supervisor's implementation

This pipeline follows the supervisor's `showershapereweighting` package
(`tree.C` + `var.h`) as closely as possible:

| Aspect | This code | Supervisor's code | Note |
|--------|-----------|-------------------|------|
| TProfile bins | `EBins` from `var.h` | Same | Matched exactly |
| Eta binning | 14 bins, [0, 2.4] | Same | Matched exactly |
| Energy function | Same `Energy(eta,phi,cells)` | Same | Same sub-window logic |
| Correction fill | `Fill(f_MC, f_PS − f_MC)` | `Fill(x, alpha)` | Identical |
| Correction apply | `Interpolate(f_MC)` | `Interpolate(MC_fraction)` | Identical |
| Cell data type | `double` (MC23e) | `float` (Run 2) | MC23e stores doubles |
| $w_{\eta\,2}$ | Relative eta positions | Absolute eta from `clusterEta` | MC23e lacks `clusterEta`; see below |
| Negative cells | Allowed (no clamping) | Same | Both allow negative cell energies |
| Conversion type | Not binned | 3 types (Inc/NCon/Con) | To be added for real data |

### Note on $w_{\eta\,2}$ calculation

The supervisor's `calcWeta()` uses absolute cell $\eta$ positions from the
`clusterEta` branch. MC23e ntuples do not have this branch, so we use
relative positions $(-1, 0, +1) \times \Delta\eta$. The width calculation is
equivalent because $w_{\eta\,2}$ is defined as a second-moment width, and the
subtraction $\langle\eta^2\rangle - \langle\eta\rangle^2$ is invariant under
translation (the absolute offset cancels).
