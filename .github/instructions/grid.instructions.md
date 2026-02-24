---
applyTo: "grid/**"
---

# Grid Submission Scripts

## Required Structure

Every grid submission script must include:

1. **Full environment setup** (do not assume anything is pre-loaded):
   ```bash
   export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
   source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet
   asetup Athena,25.0.40
   source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh"
   lsetup panda rucio
   ```

2. **Grid proxy check**:
   ```bash
   voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas
   ```

3. **Submit from `run/` directory**:
   ```bash
   cd "${WORKDIR}/run"
   ```

## GRL Verification

Before adding a new dataset, verify the run is in the GRL:
```bash
grep "<Run>RUNNUMBER</Run>" /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/<path>.xml
```

GRL paths (from NTupleMakerConfig.py `getGRL()`):
- 2022: `data22_13p6TeV/20250321/...periodAllYear...v134-pro28-09...xml` (runs 431810–440613)
- 2023: `data23_13p6TeV/20250321/...periodAllYear...v133-pro31-11...xml`
- 2024: `data24_13p6TeV/20241118/physics_25ns_data24.xml` (runs 473235–486706)

## Naming Convention

- Script: `submit_<dataset>_<campaign>.sh` (e.g., `submit_data24.sh`, `submit_mc23e_Zllg.sh`)
- Output dataset: `user.femarta.<runID>.<tag>.<version>`
- Version: `v$(date +%Y%m%d_%H%M)` for unique resubmission names
