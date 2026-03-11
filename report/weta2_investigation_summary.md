# w_eta_2 Investigation Summary

## Background

The shower shape variable $w_{\eta\,2}$ (width in eta of the shower in Layer 2) computed from NTupleMaker cell-level branches showed a discrepancy versus the stored reconstruction-level value. While $R_\eta$ and $R_\phi$ (energy ratios) showed good closure, $w_{\eta\,2}$ did not match.

This document summarises the investigation: what was tried, what worked, what failed, and possible next steps.

---

## The Variable

$$w_{\eta\,2} = \sqrt{\frac{\sum_i E_i\,\eta_i^2}{\sum_i E_i} - \left(\frac{\sum_i E_i\,|\eta_i|}{\sum_i E_i}\right)^{\!2}}$$

Computed over a **3×5 (eta × phi) window** in Layer 2 cells centred on the hottest cell.

Key difference from $R_\eta$ and $R_\phi$: those are pure energy-sum ratios (e.g. $E_{3\times7}/E_{7\times7}$) and do not depend on cell positions. $w_{\eta\,2}$ requires knowledge of **each cell's physical eta coordinate** — this is where the problems arise.

---

## What Was Investigated

### 1. Cell Ordering in NTupleMaker (-d 0)

**Finding:** When running with `-d 0` on AODs, `egammaCellDecorator` uses `CaloFillRectangularCluster` to create a 7×11 cell grid around the cluster centre. The cells are stored in **eta-major order**: index = `ie * 11 + ip`, producing a spatially-ordered rectangular array.

**However**, the earlier v1 ntuples (produced on the grid) stored cells in hash order (from `getCellLinks()` iteration), not in the proper eta×phi grid. This meant cell indices did not correspond to physical positions.

**Status:** This was identified as a key issue. New ntuples from master branch (v3 grid job) use the proper `CaloFillRectangularCluster` ordering.

### 2. Missing Cell Eta/Phi Branches

**Finding:** NTupleMaker's master branch already computed cell eta and phi positions internally (stored in `m_cells_eta` and `m_cells_phi` maps in `fillClusterCells()`), but never wrote them to the output tree. Only cell energies (`m_cells_e`) were written as branches.

**Fix:** Added 4 lines to `connectCells()` in NTupleMaker.h to create tree branches for `%ix%iClusterLr%iEta` and `%ix%iClusterLr%iPhi` (commit `044960c`). Now the cell eta/phi positions are available in the ntuple output.

### 3. ATLAS Reconstruction Geometric Correction

**Finding:** The ATLAS reconstruction applies a **geometric correction** (`egammaqweta2c::Correct()`) to the raw $w_{\eta\,2}$ value. This is a polynomial correction in $|\eta|$ that accounts for:
- Finite cell size effects
- Cell position quantisation
- Barrel/endcap geometry differences

The correction is **not** a simple scaling — it's a region-dependent polynomial transformation. The stored `weta2` in the ntuple is the **corrected** value, not the raw cell-level computation.

**Status:** The exact polynomial coefficients and correction formula were traced in the ATLAS reconstruction source code. Applying `egammaqweta2c::Correct()` in the analysis script improved agreement but did not fully close the gap.

### 4. Per-Eta-Region Analysis

**Finding:** The discrepancy between computed and stored $w_{\eta\,2}$ varies by detector region:
- **Barrel** (|η| < 0.8): Smallest discrepancy
- **Transition** (0.8 < |η| < 1.37): Moderate discrepancy
- **Crack** (1.37 < |η| < 1.52): Very few events (excluded)
- **Endcap** (1.52 < |η| < 2.47): Larger discrepancy, likely due to different cell granularity

**Status:** Documented but not fully resolved.

### 5. EGAM3 DAOD Attempt

**Finding:** Tried to use Luca's DAOD_EGAM3 files (which are smaller than AODs) to avoid grid processing of large AODs. Two approaches failed:
- `-d 1` (DAOD mode): Cell decorations (`cells_E`, `cells_eta`, etc.) do not exist in DAOD_EGAM3 → `ncells=0` for all events.
- `-d 0` (AOD mode on DAOD): `AllCalo` container not available → FPE/DIVBYZERO.

**Status:** Blocked. See [egam3_problem_report.md](egam3_problem_report.md) for full details and recommendations.

---

## What Worked

| Item | Detail |
|------|--------|
| $R_\eta$ computation | $E_{3\times7}/E_{7\times7}$ matches stored value — only depends on energy sums |
| $R_\phi$ computation | $E_{3\times3}/E_{3\times7}$ matches stored value — only depends on energy sums |
| L2 eta/phi branches | Added to NTupleMaker, now available in ntuple output |
| AOD processing with `-d 0` | Works correctly on master branch, produces complete cell data |
| Grid production v3 | Task 49065711 running on full MC23e AODs (1827 files, ~35M events, 100% on disk) |
| Closure test pipeline | Operates at the cell-energy level, $w_{\eta\,2}$ closure works because corrections propagate through the formula |

---

## What Failed / Remains Unresolved

| Item | Detail |
|------|--------|
| $w_{\eta\,2}$ from cells vs stored | Still shows discrepancy after applying `egammaqweta2c::Correct()` |
| Cell ordering in v1 ntuples | Hash-ordered cells meant 3×5 window couldn't be reliably identified from index |
| DAOD_EGAM3 | Completely blocked — no cell-level access (see separate report) |
| Exact correction coefficients | Uncertain whether the polynomial coefficients in our analysis match the Athena version used for reconstruction of these samples |

---

## Current Status

- **Grid job v3** (task 49065711) is processing full MC23e AODs on the master branch with L2 cell eta/phi branches
- Once output is available: cells will be in proper rectangular grid order with physical eta/phi positions as explicit branches
- The closure test and data-MC comparison pipelines work correctly at the cell-energy level and do not depend on resolving the $w_{\eta\,2}$ computation discrepancy

---

## Possible Next Steps / Ideas

1. **Use L2 eta/phi branches for proper window selection**
   - With explicit cell eta/phi now available, identify the hottest cell by energy, then select the 3×5 neighbours by physical proximity (eta/phi sorting) rather than relying on array index ordering.
   - This approach is robust to any cell ordering in the ntuple.

2. **Step-by-step comparison with ATLAS reconstruction**
   - Run the ATLAS reconstruction's `calcWeta()` (from `egammaMiddleShape.cxx`) on the same events and compare intermediate values at each step.
   - Check specifically: which cell is identified as the "hottest cell", and do the 3×5 window boundaries match?

3. **Verify `egammaqweta2c::Correct()` coefficients**
   - The correction polynomial coefficients may differ between Athena versions. Verify that the coefficients used in our analysis match those in Athena 25.0.40 (the version used for reconstruction of these MC samples).
   - The relevant source: `Reconstruction/egamma/egammaUtils/src/egammaqweta2c.cxx`.

4. **Supercluster vs sliding-window cluster**
   - ATLAS reconstruction may compute $w_{\eta\,2}$ on the supercluster (topo-cluster based) rather than the sliding-window cluster. Check if NTupleMaker's `egammaCellDecorator` uses the same cluster type.
   - The cluster type affects which cells are in the 3×5 window.

5. **Negative energy cells**
   - Check whether ATLAS reconstruction excludes negative-energy cells from the $w_{\eta\,2}$ sum. NTupleMaker currently includes all cells.

6. **For DAOD_EGAM3**
   - Report to Luca/supervisor that the derivation needs a cell decorator algorithm to write per-object cell vectors.
   - Alternatively, continue using AODs with `-d 0` — the 35M events from the v3 grid job provide ample statistics.
