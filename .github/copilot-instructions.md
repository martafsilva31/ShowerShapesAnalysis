# ShowerShapesAnalysis — Copilot Instructions

> **Scope**: This file applies to `QT/ShowerShapes/ShowerShapesAnalysis/` and its contents.
> It supplements (does NOT repeat) the global instructions at `../../.github/copilot-instructions.md`.
> See the global file for: multi-agent workflow, self-maintenance discipline, operating principles, task tracking.

---

## 1. Project Overview

ATLAS photon shower shape reweighting analysis. Computes and validates calorimeter shower shape variables (R_eta, R_phi, w_eta_2) from Layer 2 cell energies, compares fudged vs unfudged values, and applies cell-based reweighting corrections.

**Physics goal**: Derive data-driven corrections for photon shower shapes used in photon ID, by comparing data and MC at the cell level in the second electromagnetic calorimeter layer.

---

## 2. Architecture

```
ShowerShapesAnalysis/          # This repo (git → github.com/martafsilva31/ShowerShapesAnalysis)
├── .github/                   # Copilot instructions (this file + scoped instructions)
│   ├── copilot-instructions.md
│   └── instructions/
│       ├── scripts.instructions.md   # Scoped to scripts/**
│       └── grid.instructions.md      # Scoped to grid/**
├── scripts/                   # Analysis scripts
│   ├── utils.py               # Shared I/O and cell-energy utilities
│   ├── compute_reta.py / compute_rphi.py / compute_weta_2.py   # Shower shape computations
│   ├── plot_reta.C / plot_rphi.C / plot_weta_2.C               # ROOT plot macros
│   ├── layer2_reweighting/          # ★ ACTIVE: Layer 2 cell-energy reweighting pipeline
│   │   ├── config.h           # Shared config: cuts, branches, geometry, formulas, pT bins
│   │   ├── fill_histograms.C  # Two-pass pipeline: accumulate → correct → fill
│   │   ├── plot_shower_shapes.C # Shower shape plots: per-eta + per-pT PDFs
│   │   ├── plot_cell_profiles.C # 7×11 cell heatmaps + correction vector plots
│   │   ├── plot_mu_scaling.C  # Pileup study: χ²/ndf and M1 shift vs ⟨μ⟩
│   │   ├── extract_chi2.C     # Chi-squared tables for a single variant
│   │   ├── extract_comparison_final.C  # Cross-variant chi2 comparison table
│   │   ├── run_layer2_final.sh  # Full pipeline driver (6 phases including compendium)
│   │   ├── make_compendiums.py  # LaTeX compendium generator for all variants/scenarios
│   │   └── README.md
│   ├── layer1_reweighting/          # ★ ACTIVE: Layer 1 cell-energy reweighting
│   │   ├── config.h / fill_histograms.C / plot_shower_shapes.C / plot_cell_profiles.C
│   │   ├── extract_chi2.C
│   │   ├── run_layer1_final.sh
│   │   └── README.md
│   ├── data_mc_normalisation_mismatches/  # ★ ACTIVE: Data–MC kinematic normalisation study
│   │   ├── plot_kinematics.C / plot_kinematics_abs.C
│   │   ├── run_kinematics.sh / run_kinematics_phase2.sh
│   │   ├── make_kinematics_compendium.py
│   │   └── README.md
│   └── old/                         # Archived: closure_test/, cell_collision_study/
├── grid/                      # Grid submission scripts (pathena)
│   ├── samples/               # Dataset sample lists
│   ├── download_ntuples.sh    # Download + merge ntuples from grid
│   ├── submit_data24.sh       # Data 2024
│   ├── submit_data24_all.sh / submit_data24_egam3.sh / submit_data24_egam4.sh
│   ├── submit_mc23e_DAOD.sh / submit_mc23e_Zllg.sh
│   ├── resubmit_477048_egam3.sh
│   └── old/                   # Archived: submit_data22.sh (Run 2, GRL mismatch — superseded)
├── configs/                   # Analysis configuration files (gitignored)
├── ntuples/                   # gitignored — ROOT ntuples
│   ├── data24/                # Data 2024
│   └── mc23e/                 # mc23e_700770_Zeeg.root (13 GB), mc23e_700771_Zmumug.root (23 GB)
├── output/                    # Mostly gitignored — computed histograms + plots
│   └── Layer_2/               # Current output for Layer 2 cell reweighting
│       ├── eta_loose/ eta_tight/ eta_pt_loose/ eta_pt_tight/ eta_mu_loose/
│       │   └── {channel}/{scenario}/
│       │       ├── histograms.root
│       │       ├── result_compendium_{variant}_{channel}_{scenario}.pdf  ← tracked in git
│       │       └── plots/
├── report/                    # LaTeX reports and documentation
│   ├── layer2/                    # cell_reweighting_report.tex
│   ├── data_mc_normalisation_mismatches/  # kinematics_report.tex
│   ├── old/                       # Archived: closure test, theory, egam3, weta2, comparison_report/
│   └── refs_theory.bib
└── setup.sh                   # Environment: setupATLAS + Athena 25.0.40 + NTupleMaker build
```

### Sibling Repos (same parent `QT/ShowerShapes/`)

| Repo | Path | Remote | Purpose |
|------|------|--------|---------|
| NTupleMaker | `../NTupleMaker_workspace/source/NTupleMaker/` | `gitlab.cern.ch/femarta/cellntuplemaker` | Athena C++ algorithm producing cell-level ntuples |
| showershapereweighting | `../showershapereweighting/` | `gitlab.cern.ch/mobelfki/showershapereweighting` | Supervisor's Run 2 correction method (CalCoef.C, NTUP.C) |

---

## 3. Technical Context

### Reference Note

**ATL-COM-PHYS-2021-640** (located at `QT/ShowerShapes/ATL-COM-PHYS-2021-640.pdf`) is the foundational document for this analysis. It describes two methods for correcting data-MC shower shape disagreements using the full Run 2 dataset:

1. **Fudge factors** (Section 4): Shift (and optionally stretch) MC shower shape distributions to match data. A χ² scan compares smoothed (KDE) PDFs of data and MC in bins of pT, |η|, and conversion type. Shift-only: `SSnew = SSold + shift`. Shift+stretch: `SSnew = stretch × (SSold − stretch_point) + shift + stretch_point`.
2. **Cell-energy reweighting** (Section 5): Correct MC cell energies in the 2nd EM calorimeter layer, so that all shower shapes built from those cells are improved simultaneously — a lower-level correction than fudge factors.

### Cell-Energy Reweighting Method (from note Section 5.2)

The **new cell reweighting method** (proposed in the note) corrects both the mean and variance of normalized cell energy distributions by applying shifts and stretches at each cell:

1. For each cell *i* in a 7×11 (η×φ) cluster, compute normalized energy: `e_i = E_i / E_cluster`
2. Compute mean (`ē`) and RMS of normalized energy distributions for data and MC.
3. Reweighted normalized cell energy:
   ```
   e_i^{MC-RW} = (RMS_data_i / RMS_MC_i) × (e_i^MC − ē_i^MC) + ē_i^data
   ```
4. Reweighted cell energy (preserving cluster energy):
   ```
   E_i^{MC-RW} = (RMS_data_i / RMS_MC_i) × E_i^MC + (ē_i^data − (RMS_data_i / RMS_MC_i) × ē_i^MC) × E^MC
   ```
5. Final rescaling: `E_i^{MC-RW} *= Σ E_i^MC / Σ E_i^{MC-RW}` to guarantee constant cluster energy.

This yields two correction matrices per bin: **shift** and **stretch** (one value per cell in the 7×11 cluster).

**Previous method (Section 5.1)**: Only corrected the mean (shift-only): `E_i^{MC-RW} = E_i^MC + Δ_i × E^MC`, where `Δ_i = ē_i^data − ē_i^MC`. This method (derived for electrons) did not work well for photons.

### Binning Strategies (from note Section 5.3)

Corrections are parameterized in several strategies:
- **Strategy 1**: cluster energy (E_cluster)
- **Strategy 2**: pT and |η|
- **Strategy 3**: pT and λ (radial distance from topo-cluster barycenter to central cell center: `λ = sqrt(Δη² + Δφ²)`)
- **Strategy 4**: pT and 4 quadrants of central cell (sign of Δη and Δφ)
- **Strategy 5**: pT and combined λ + quadrants (centered region λ < 6.25×10⁻³, plus 4 quadrants for off-center hits)

Standard binning:
- pT: [10, 15, 20, 25, 30, 40, 1000] GeV
- |η|: [0, 0.6, 0.8, 1.37, 1.52, 1.81, 2.01, 2.37]
- λ: [0, 3.125, 6.250, 9.750, 17.677] × 10⁻³
- Corrections derived separately for unconverted and converted photons.

### Key Results (from note Section 5.4 & Section 6)

- Cell reweighting achieves excellent agreement between data and MC for R_eta, R_phi, and w_eta_2 — often outperforming fudge factors (TUNE22).
- Advantage: cell-level corrections automatically improve ALL shower shapes, unlike fudge factors which correct each variable independently.
- Limitation (at time of note): only Layer 2 corrections derived; extension to all layers planned.
- Ideal approach: cell reweighting for bulk corrections, fudge factors for residual differences.

### Event Selection for Cell Reweighting (from note Section 3.5)

- Only "healthy clusters": 77 cells (no missing cells), central cell must have highest energy ("hottest cell").
- No photon ID or FixedCutLoose applied for cell reweighting studies.
- Angular separation photon–lepton: ΔR > 0.4.
- Pileup reweighting applied; MC reweighted to data luminosity.

### Shower Shape Variables

| Variable | Formula | Cell Window | Branch |
|----------|---------|-------------|--------|
| R_eta | E237 / E277 | 3×7 in 7×7, 7×7 total | `photon.reta`, `photon.unfudged_reta` |
| R_phi | E233 / E237 | 3×3 in 3×7 | `photon.rphi`, `photon.unfudged_rphi` |
| w_eta_2 | Energy-weighted η width | 3×5 (η×φ) | `photon.weta2`, `photon.unfudged_weta2` |

Additional shower shapes from the note (not currently computed in scripts but used for fudge factors):
- R_had / R_had1: hadronic leakage ratios
- w_eta_1 (w_1): lateral width in Layer 1 strips
- w_s_tot: total lateral width in a window of Δη ≈ 0.0625
- f_side: energy fraction outside 3 central strip cells within 7 cells
- ΔE: energy difference between 1st and 2nd maximum in strips
- E_ratio: (E_max1 − E_max2) / (E_max1 + E_max2) in strips
- f_1: fraction of energy in Layer 1

### Fudge Factor Method Details (from note Section 4)

- PDFs smoothed using adaptive KDE (TMVA) with tunable fine factors per variable.
- Binning: 6 pT bins × 8 |η| bins × conversion type.
- Shift-only: MC PDF shifted ±50 bins (out of 500) → χ² scan → parabolic fit at minimum → shift value.
- Shift+stretch: 2D scan (shift × stretch plane, stretch 0.7–1.7 in 150 steps) → minimum bin gives both corrections.
- Statistical uncertainties: shift from χ²_min + 1; shift+stretch from χ²_min + 2.3 contour (1σ ellipse).

### Cluster Cell Arrangement

- The 7×11 cluster is stored as `photon.7x11ClusterLr2E` — a `vector<double>` of 77 values in **row-major** order (7 eta × 11 phi).
- "Fudged" = after ATLAS correction factors are applied; "unfudged" = raw from calorimeter.
- Cell numbering (from note Figure 5): cell index = phi + 11 × eta (phi: 0–10, eta: 0–6). Central cell (hottest) is cell 38 in 0-based indexing (eta=3, phi=5). NOTE: the note uses 1-based indexing and calls it cell 39.
- Layer 2 cell size: Δη = 0.025, Δφ = 0.0245.

### NTupleMaker Key Code Paths

- `NTupleMakerConfig.py`: `GRLCfg()`, `getGRL(year)`, `getZllgTriggers(year)`, `NTupleMakerCfg(flags)`
- `NTupleMaker.cxx`: GRL filtering at line 155 (`m_GRLTool->passRunLB()`), cutflow bins
- `jobConfig.py`: argparse CLI with `-y year`, `-a algo`, `-d isDAODCells`, `--normalise`

### Known Issues

- **Data22 has 0 events**: Runs 428648–429027 (periods A–E) are not in the GRL which starts at run 431810 (period F). Need data from periods F/H/J.
- **GRL paths are hardcoded** in `NTupleMakerConfig.py` via `getGRL()` function, pointing to `/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/`.

### Pipeline Histogram Naming Convention

All histograms in `histograms.root` use the following naming pattern:
- Integrated: `h_{var}_{tag}` (e.g. `h_reta_data_stored`, `h_reta_mc_fudged`)
- Per-eta bin: `h_{var}_{tag}_eta{NN}` (e.g. `h_reta_data_eta00`, `h_reta_mc_M1_eta07`)
- Per-(eta,pT) bin (eta_pt mode only): `h_{var}_{tag}_eta{NN}_pt{PP}` (e.g. `h_reta_data_eta00_pt02`)
- Cell profiles (in `cell_profiles/` subdir): `h_frac_{type}_eta{NN}` or `h_frac_{type}_eta{NN}_pt{PP}`
- Corrections (in `corrections/` subdir): `h_delta_eta{NN}` or `h_delta_eta{NN}_pt{PP}`

Variable names: `reta`, `rphi`, `weta2`. Tags: `data`, `data_stored`, `mc`, `mc_fudged`, `mc_unfudged`, `mc_M1`, `mc_M2`.

### Pipeline Variants

| Variant | Binning | Isolation | Correction dimensions |
|---------|---------|-----------|----------------------|
| `eta_loose` | 14 eta bins | Loose (default) | `[14][1][77]` |
| `eta_tight` | 14 eta bins | Tight | `[14][1][77]` |
| `eta_pt_loose` | 14 eta × 6 pT bins | Loose | `[14][6][77]` |
| `eta_pt_tight` | 14 eta × 6 pT bins | Tight | `[14][6][77]` |
| `eta_mu_loose` | 14 eta × 4 ⟨μ⟩ bins | Loose | `[14][4][77]` |

pT bins (GeV): [10, 15, 20, 25, 30, 40, 1000]. Defined in `config.h` as `kPtLimits`.
⟨μ⟩ bins: [0, 40, 55, 70, 120]. Defined in `config.h` as `kMuLimits`.

`plot_shower_shapes.C(channel, scenario, baseDir, binning="eta", isolation="loose")` produces:
- `rew_{reta,rphi,weta2}.pdf` — per-eta shower shape comparison (Data, MC, M1, M2)
- `rew_integrated.pdf` — integrated (all eta combined)
- `computed_vs_stored.pdf`, `computed_vs_stored_eta.pdf` — cell-computed vs branch values
- `fudge_factors.pdf`, `fudge_factors_eta.pdf` — fudge factor comparison
- When `binning="eta_pt"`: `rew_{var}_pt{PP}.pdf` — per-pT shower shape comparison
- When `binning="eta_mu"`: `rew_{var}_mu{PP}.pdf` — per-⟨μ⟩ shower shape comparison

`plot_cell_profiles.C(channel, scenario, baseDir, binning="eta", isolation="loose")` produces:
- `cell_{data,mc,mc_m1,mc_m2}.pdf` — 7×11 cell heatmaps per eta
- `cell_{shift,stretch}.pdf` — correction vector maps per eta
- When `binning="eta_pt"`: `cell_{type}_pt{PP}.pdf` — per-pT cell heatmaps
- When `binning="eta_mu"`: `cell_{type}_mu{PP}.pdf` — per-⟨μ⟩ cell heatmaps

`plot_mu_scaling.C(baseDir, channel, scenario)` — for `eta_mu` variants only — produces:
- `plots/mu_chi2_scaling.pdf` — χ²/ndf vs ⟨μ⟩ bin for MC/M1/M2, one page per |η| bin
- `plots/mu_m1_shift_scaling.pdf` — mean |M1 shift| per cell vs ⟨μ⟩ bin, one page per |η| bin

Output base directory: `output/Layer_2/{variant}/{channel}/{scenario}/`
where `variant` is one of `eta_loose`, `eta_tight`, `eta_pt_loose`, `eta_pt_tight`, `eta_mu_loose`.
Scripts use relative path `../../output/Layer_2/...` from `scripts/layer2_reweighting/`.

---

## 4. Conventions for This Repo

### Scripts (`scripts/`)
- CLI-based: accept file paths as arguments, no hardcoded paths.
- Naming: `<verb>_<observable>.py` (e.g., `compute_weta2.py`, `plot_reta_rphi.C`).
- Include a docstring or header comment with description, inputs, outputs.
- Use ROOT's `TFile`, `TTree`, `TH1` — no custom frameworks.

### Grid (`grid/`)
- Every script must include full env setup (setupATLAS, asetup, lsetup, build source, voms-proxy check).
- Must `cd` to `NTupleMaker_workspace/run/` before calling `pathena`.
- Verify runs are in GRL before creating new submission scripts.
- Dataset naming: `user.femarta.<runID>.<tag>.<version>`.

### Git
- Branch from `master` for features.
- Conventional commits: `feat:`, `fix:`, `docs:`, `chore:`.
- Never commit: `.root`, `ntuples/`, `plots/`, `*.log`, `*.npz`, `*.pdf`, `*.png`.
- Propose exact `git add`, `git commit -m "..."`, `git push` commands.
- **Commit and push after each completed phase** (fill, plot, report, docs). Never accumulate more than one phase in a single commit.

---

## 5. Documentation Discipline (Specific to This Repo)

In addition to the global self-maintenance rules, for THIS repo:

| When This Happens | Update These |
|-------------------|-------------|
| New script added to `scripts/` | README.md § Directory Structure + this file § Architecture |
| New grid script added | README.md § Grid Submission + `grid.instructions.md` |
| New ntuple campaign downloaded | README.md § Known Issues (if relevant), this file § Architecture |
| Active study folder reorganised | README.md § Active Studies + this file § Architecture |
| New known issue found | README.md § Known Issues + this file § Known Issues |
| Setup procedure changes (new Athena version, new tool) | `setup.sh` + this file § Architecture + global instructions § Software Stack |

**Rule**: A task is NOT complete until the README and copilot instructions are consistent with the code.
