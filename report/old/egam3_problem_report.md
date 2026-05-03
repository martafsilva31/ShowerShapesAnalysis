# DAOD_EGAM3 Problem Report — NTupleMaker Cell-Level Access

## Summary

DAOD_EGAM3 files (from Luca's derivation) **cannot be used** with NTupleMaker to extract individual calorimeter cell energies and positions. Neither of NTupleMaker's two code paths (`-d 0` or `-d 1`) works on DAOD_EGAM3. The current workaround is to run on AOD files with `-d 0`.

---

## Background

NTupleMaker extracts photon/electron calorimeter cell information to compute shower shape variables (R_eta, R_phi, w_eta_2) from individual cell energies. It has two cell-reading code paths controlled by the `-d` flag:

| Flag | Mode | Input type | How cells are read |
|------|------|------------|--------------------|
| `-d 0` | AOD mode | AOD | `egammaCellDecorator` tool reads the `AllCalo` cell container |
| `-d 1` | DAOD mode | DAOD | Reads pre-existing decorations on the egamma object |

---

## Problem: `-d 1` on DAOD_EGAM3

The DAOD code path (`-d 1`) expects the following decorations on each photon/electron object:

| Decoration | Type | Purpose |
|------------|------|---------|
| `cells_E` | `vector<float>` | Cell energies |
| `cells_eta` | `vector<float>` | Cell eta coordinates |
| `cells_phi` | `vector<float>` | Cell phi coordinates |
| `cells_layer` | `vector<int>` | Cell calorimeter layer index (0–4) |
| `ncells` | `int` | Number of cells |

**These decorations do not exist in DAOD_EGAM3.** When NTupleMaker tries to read them, all vectors are empty and `ncells = 0` for every event. No cell data is available.

**Root cause:** Luca's EGAM3 derivation configuration does not include a cell decorator algorithm that writes these per-object cell vectors. The DAOD only contains cluster-level quantities (e.g. `ethad`, `Reta`, `weta2`), not the underlying cell-level data.

---

## Problem: `-d 0` on DAOD_EGAM3

The AOD code path (`-d 0`) uses the `egammaCellDecorator` tool, which:

1. Reads the full detector cell container (`AllCalo` / `CaloCellContainer`)
2. Finds cells around the cluster centre by looping over the calorimeter neighbourhood
3. Writes per-object decorations: `cell_%ix%iClusterLr%iE`, `cell_%ix%iClusterLr%iEta`, etc.

**This fails on DAOD_EGAM3** because DAODs do not include the full `AllCalo` container. DAODs only store a thinned subset of cells (within a narrow cone around the photon), and the full cell topology required by `egammaCellDecorator` is not available. Running `-d 0` on the DAOD produces floating-point exceptions (FPE/DIVBYZERO) and unreliable output.

---

## What Would Fix It

To make DAOD_EGAM3 usable for cell-level analysis, the derivation configuration needs to include a cell decorator algorithm that writes the required per-object decorations. Specifically:

1. **Add a cell decoration step** to the EGAM3 derivation framework configuration. This could be:
   - A `CellDecoratorAlg` (or equivalent Athena algorithm) that reads the thinned cell container and writes `cells_E`, `cells_eta`, `cells_phi`, `cells_layer`, `ncells` vectors as decorations on each photon/electron.
   - Or adapt `egammaCellDecorator` to work with the DAOD's thinned cell container instead of requiring the full `AllCalo`.

2. **Ensure the thinned cell container retains enough cells.** The current EGAM3 thinning may discard cells outside a narrow cone. For shower shape computation, we need at least the full 7×11 cell window (in Layer 2 eta × phi units) around the cluster centre.

3. **Alternative**: Stick with AOD files and `-d 0`. AODs are larger (~7 GB/file vs ~1 GB for DAODs) but contain the full `AllCalo` container and work correctly with NTupleMaker.

---

## Current Workaround

We are running on **AOD files** with `-d 0` (AOD mode). This works correctly:

- Grid job task **49065711** is processing the full MC23e eegamma AOD dataset
- Input: `mc23_13p6TeV.700770.Sh_2214_eegamma.merge.AOD.e8514_e8586_s4369_s4370_r16083_r15970` (tid42517418, 100% on disk at FZK-LCG2)
- 1827 files, ~35M events, 1 file per job
- NTupleMaker on `master` branch with L2 cell eta/phi branches added (commit `044960c`)

---

## Test Evidence

| Test | Command | Result |
|------|---------|--------|
| DAOD_EGAM3 + `-d 1` | `NTupleMaker.jobConfig -d 1` on DAOD_EGAM3 file | `ncells=0` for all events, empty cell vectors |
| DAOD_EGAM3 + `-d 0` | `NTupleMaker.jobConfig -d 0` on DAOD_EGAM3 file | FPE/DIVBYZERO warnings, unreliable output |
| AOD + `-d 0` | `NTupleMaker.jobConfig -d 0` on AOD file | **Works correctly**: 4112 entries from 30k events, all branches populated |

---

## Relevant Code Locations (NTupleMaker)

- `egammaCellDecorator` initialization guard: `NTupleMaker.cxx`, line 64 — `if(not m_isDAODCells) CHECK(m_egammaCellDecorator.retrieve())`
- Cell decoration call: `NTupleMaker.cxx`, line 324 — `if(not m_isDAODCells) CHECK(m_egammaCellDecorator->decorate((*photon)))`
- DAOD decoration reading: `NTupleMaker.cxx`, `fillClusterCells()` function (~line 543) — reads `cells_E`, `cells_eta`, `cells_phi`, `cells_layer`, `ncells`
- `egammaCellDecorator.cxx` — reads from `AllCalo` container (`CaloCellContainer`), produces `cell_%ix%iClusterLr%iE/Eta/Phi` decorations
