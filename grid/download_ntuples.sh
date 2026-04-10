#!/bin/bash
###############################################################################
# download_ntuples.sh
#
# Downloads all full-production ntuples (EGAM3 + EGAM4 + MC23e) from
# completed grid jobs, then merges them into one ROOT file per category
# and one combined file.
#
# Output directory: /dcache/atlas/mfernand/qt_ntuples/data24/
#
# Output structure:
#   raw/egam3/          ← rucio downloads (one subdir per run)
#   raw/egam4/
#   raw/mc_eegamma/
#   raw/mc_mumugamma/
#   egam3.root          ← merged data Zeeg (195 runs)
#   egam4.root          ← merged data Zmumug (195 runs)
#   mc_eegamma.root     ← merged MC Zeeg (700770, v1r2)
#   mc_mumugamma.root   ← merged MC Zmumug (700771, v1)
#   all_merged.root     ← all four categories combined
#
# Usage:
#   ./download_ntuples.sh                  # download + merge all
#   ./download_ntuples.sh --skip-download  # merge only (files already local)
#   ./download_ntuples.sh --merge-only     # same as --skip-download
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTDIR="/dcache/atlas/mfernand/qt_ntuples/data24"
RAWDIR="${OUTDIR}/raw"
RUCIO_ACCOUNT="${RUCIO_ACCOUNT:-femarta}"
NDOWNLOADERS=16    # parallel streams per rucio download call

# Colours
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

SKIP_DOWNLOAD=false
for arg in "$@"; do
    [[ "$arg" == "--skip-download" || "$arg" == "--merge-only" ]] && SKIP_DOWNLOAD=true
done

echo -e "${YELLOW}============================================================${NC}"
echo -e "${YELLOW}  Full-Production Ntuple Download  (EGAM3 + EGAM4 + MC23e)${NC}"
echo -e "${YELLOW}============================================================${NC}"
echo "  Output: ${OUTDIR}"
echo ""

# ---------------------------------------------------------------------------
# ATLAS environment setup (matches working submit scripts)
# ---------------------------------------------------------------------------
echo -e "${GREEN}[0/5] Setting up ATLAS environment...${NC}"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
# ATLAS setup scripts use unbound variables — suspend strict mode
set +eu
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet || true
# NOTE: Do NOT use 'asetup Athena' here — it conflicts with rucio and causes
# silent download failures (empty directories). Only rucio + ROOT are needed.
export RUCIO_ACCOUNT=femarta
lsetup rucio || true
lsetup "root 6.28.04-x86_64-el9-gcc13-opt" --quiet || true
set -eu
echo "  rucio account: ${RUCIO_ACCOUNT}"
echo ""

# Grid proxy
if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
    export X509_VOMS_DIR=/cvmfs/grid.cern.ch/etc/grid-security/vomsdir
    if ! voms-proxy-info --exists &>/dev/null; then
        echo -e "${YELLOW}  No valid proxy found.${NC}"
        echo -e "${YELLOW}  Enter your GRID CERTIFICATE PASSWORD when prompted:${NC}"
        echo ""
        voms-proxy-init -voms atlas
        echo ""
    else
        echo "  Grid proxy: $(voms-proxy-info --timeleft)s remaining"
    fi
    echo ""
fi

# Sanity checks
if [[ "${SKIP_DOWNLOAD}" == false ]] && ! command -v rucio &>/dev/null; then
    echo -e "${RED}ERROR: rucio not in PATH after lsetup.${NC}"; exit 1
fi
if ! command -v hadd &>/dev/null; then
    echo -e "${RED}ERROR: hadd not in PATH after ROOT lsetup.${NC}"; exit 1
fi

# Create raw subdirectories
mkdir -p "${RAWDIR}/egam3" "${RAWDIR}/egam4" "${RAWDIR}/mc_eegamma" "${RAWDIR}/mc_mumugamma"

# ---------------------------------------------------------------------------
# Step 1: Build dataset lists from BigPandaMon
# ---------------------------------------------------------------------------
echo -e "${GREEN}[1/5] Querying BigPandaMon for best output datasets...${NC}"

TMPDIR_LISTS="$(mktemp -d)"
trap 'rm -rf "${TMPDIR_LISTS}"' EXIT

python3 - "${TMPDIR_LISTS}" << 'PYEOF'
import sys, urllib.request, json, ssl, re

outdir = sys.argv[1]
ctx = ssl.create_default_context()

SCRIPT_BASE = "/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace/source/NTupleMaker/script"

def get_runs(filename):
    runs = []
    with open(f"{SCRIPT_BASE}/{filename}") as f:
        for line in f:
            m = re.search(r'# Run (\d+)', line)
            if m:
                runs.append(m.group(1))
    return runs

egam3_runs = get_runs("data24_13p6TeV.DAOD_EGAM3.25ns.p61xx.txt")
egam4_runs = get_runs("data24_13p6TeV.DAOD_EGAM4.25ns.p61xx.txt")

url = 'https://bigpanda.cern.ch/tasks/?taskname=user.femarta.*&json&limit=500&days=30'
req = urllib.request.Request(url, headers={'Accept': 'application/json'})
resp = urllib.request.urlopen(req, context=ctx, timeout=60)
tasks = json.loads(resp.read())

def extract_run(taskname):
    m = re.search(r'user\.femarta\.(\d+)\.', taskname)
    return m.group(1) if m else None

categories = {
    'egam3':       'Zeeg_data24_egam3',
    'egam4':       'Zmumug_data24_egam4',
    'mc_eegamma':  'Zeeg_mc23e_DAOD',
    'mc_mumugamma':'Zmumug_mc23e_DAOD',
}
by_cat = {cat: {} for cat in categories}

for t in tasks:
    name = t.get('taskname', '').rstrip('/')
    for cat, pattern in categories.items():
        if pattern in name:
            run = extract_run(name)
            if run:
                by_cat[cat].setdefault(run, []).append({
                    'tid': t.get('jeditaskid'),
                    'taskname': name,
                    'status': t.get('status'),
                })

priority = {'done':0,'finished':1,'running':2,'scouting':3,
            'exhausted':5,'failed':6,'broken':7,'aborted':8,'tobroken':9}

def best(lst):
    return sorted(lst, key=lambda x: (priority.get(x['status'], 99), -x['tid']))[0]

datasets = {'egam3': [], 'egam4': [], 'mc_eegamma': [], 'mc_mumugamma': []}

for run in egam3_runs:
    padded = f"00{run}"
    tl = by_cat['egam3'].get(padded, []) or by_cat['egam3'].get(run, [])
    if tl:
        bt = best(tl)
        if bt['status'] in ('done', 'finished'):
            ver = bt['taskname'].split('.')[-1]
            datasets['egam3'].append(f"user.femarta.{padded}.Zeeg_data24_egam3.{ver}_EXT0")

for run in egam4_runs:
    padded = f"00{run}"
    tl = by_cat['egam4'].get(padded, []) or by_cat['egam4'].get(run, [])
    if tl:
        bt = best(tl)
        if bt['status'] in ('done', 'finished'):
            ver = bt['taskname'].split('.')[-1]
            datasets['egam4'].append(f"user.femarta.{padded}.Zmumug_data24_egam4.{ver}_EXT0")

for cat, dsid, prefix in [('mc_eegamma','700770','Zeeg_mc23e_DAOD'),
                           ('mc_mumugamma','700771','Zmumug_mc23e_DAOD')]:
    tl = by_cat[cat].get(dsid, [])
    if tl:
        bt = best(tl)
        if bt['status'] in ('done', 'finished'):
            ver = bt['taskname'].split('.')[-1]
            datasets[cat].append(f"user.femarta.{dsid}.{prefix}.{ver}_EXT0")

for cat, ds_list in datasets.items():
    listfile = f"{outdir}/datasets_{cat}.txt"
    with open(listfile, 'w') as f:
        f.write('\n'.join(ds_list) + '\n')
    print(f"  {cat}: {len(ds_list)} datasets -> {listfile}")

PYEOF

echo ""

# ---------------------------------------------------------------------------
# Step 2: Download with rucio (batched to avoid API overload)
# ---------------------------------------------------------------------------
BATCH_SIZE=40   # datasets per rucio download call (tested: 20 OK, 195 too many)

download_category() {
    local cat="$1"
    local rawsubdir="${RAWDIR}/${cat}"
    local listfile="${TMPDIR_LISTS}/datasets_${cat}.txt"
    local n
    n=$(wc -l < "${listfile}" 2>/dev/null || echo 0)
    n=${n//[[:space:]]/}

    if [[ "$n" -eq 0 ]]; then
        echo -e "${RED}  [${cat}] No datasets found — skipping.${NC}"
        return
    fi

    echo -e "${GREEN}  Downloading ${cat} (${n} datasets, batch=${BATCH_SIZE}) -> ${rawsubdir}/${NC}"

    local batch=() bnum=1 total_batches=$(( (n + BATCH_SIZE - 1) / BATCH_SIZE ))
    while IFS= read -r ds; do
        [[ -z "$ds" ]] && continue
        batch+=("$ds")
        if [[ "${#batch[@]}" -ge "${BATCH_SIZE}" ]]; then
            echo -e "    batch ${bnum}/${total_batches} (${#batch[@]} datasets)..."
            rucio download \
                --ndownloader "${NDOWNLOADERS}" \
                --dir "${rawsubdir}" \
                "${batch[@]}" || echo -e "${RED}    WARNING: batch ${bnum} had errors${NC}"
            batch=()
            (( bnum++ ))
        fi
    done < "${listfile}"
    # remaining datasets
    if [[ "${#batch[@]}" -gt 0 ]]; then
        echo -e "    batch ${bnum}/${total_batches} (${#batch[@]} datasets)..."
        rucio download \
            --ndownloader "${NDOWNLOADERS}" \
            --dir "${rawsubdir}" \
            "${batch[@]}" || echo -e "${RED}    WARNING: batch ${bnum} had errors${NC}"
    fi
    echo -e "${GREEN}  [${cat}] done.${NC}"
}

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    echo -e "${GREEN}[2/5] Downloading ntuples (4 categories in parallel)...${NC}"
    LOGDIR="${OUTDIR}/logs"
    mkdir -p "${LOGDIR}"
    download_category egam3       > "${LOGDIR}/egam3.log"       2>&1 &
    PID_EGAM3=$!
    download_category egam4       > "${LOGDIR}/egam4.log"       2>&1 &
    PID_EGAM4=$!
    download_category mc_eegamma  > "${LOGDIR}/mc_eegamma.log"  2>&1 &
    PID_MCEEG=$!
    download_category mc_mumugamma > "${LOGDIR}/mc_mumugamma.log" 2>&1 &
    PID_MCMUMU=$!
    echo "  PIDs: egam3=${PID_EGAM3} egam4=${PID_EGAM4} mc_eeg=${PID_MCEEG} mc_mumu=${PID_MCMUMU}"
    echo "  Per-category logs in ${LOGDIR}/"
    echo "  Waiting for all downloads to finish..."
    FAIL=0
    wait ${PID_EGAM3}  || { echo -e "${RED}  egam3 download had errors (see ${LOGDIR}/egam3.log)${NC}"; FAIL=1; }
    wait ${PID_EGAM4}  || { echo -e "${RED}  egam4 download had errors (see ${LOGDIR}/egam4.log)${NC}"; FAIL=1; }
    wait ${PID_MCEEG}  || { echo -e "${RED}  mc_eegamma download had errors (see ${LOGDIR}/mc_eegamma.log)${NC}"; FAIL=1; }
    wait ${PID_MCMUMU} || { echo -e "${RED}  mc_mumugamma download had errors (see ${LOGDIR}/mc_mumugamma.log)${NC}"; FAIL=1; }
    echo -e "${GREEN}  All downloads finished.${NC}"
    echo ""
else
    echo -e "${YELLOW}[2/5] Skipping download (--skip-download / --merge-only)${NC}"
    echo ""
fi

# ---------------------------------------------------------------------------
# Helper: merge ROOT files in a raw subdirectory
# ---------------------------------------------------------------------------
merge_category() {
    local cat="$1"
    local outfile="${OUTDIR}/${cat}.root"
    local rawsubdir="${RAWDIR}/${cat}"

    mapfile -d '' files < <(find "${rawsubdir}" -name "*.root" -type f -print0 2>/dev/null | sort -z)
    local n="${#files[@]}"

    if [[ "$n" -eq 0 ]]; then
        echo -e "${RED}  [${cat}] No ROOT files found in ${rawsubdir} — skipping merge.${NC}"
        return 1
    fi

    echo -e "${GREEN}  [${cat}] Merging ${n} file(s) -> ${outfile}${NC}"
    if [[ "$n" -eq 1 ]]; then
        cp "${files[0]}" "${outfile}"
    else
        hadd -f "${outfile}" "${files[@]}"
    fi
    echo "         Size: $(du -h "${outfile}" | cut -f1)"
    return 0
}

# ---------------------------------------------------------------------------
# Steps 3–4: Merge per category
# ---------------------------------------------------------------------------
echo -e "${GREEN}[3/5] Merging per category...${NC}"
merge_category egam3       && EGAM3_OK=true   || EGAM3_OK=false
merge_category egam4       && EGAM4_OK=true   || EGAM4_OK=false
merge_category mc_eegamma  && MCEEG_OK=true   || MCEEG_OK=false
merge_category mc_mumugamma && MCMUMU_OK=true  || MCMUMU_OK=false
echo ""

# ---------------------------------------------------------------------------
# Step 5: Combined file (all four categories)
# ---------------------------------------------------------------------------
echo -e "${GREEN}[4/5] Building combined all_merged.root...${NC}"
COMBINED_INPUTS=()
[[ "$EGAM3_OK"  == true ]] && COMBINED_INPUTS+=("${OUTDIR}/egam3.root")
[[ "$EGAM4_OK"  == true ]] && COMBINED_INPUTS+=("${OUTDIR}/egam4.root")
[[ "$MCEEG_OK"  == true ]] && COMBINED_INPUTS+=("${OUTDIR}/mc_eegamma.root")
[[ "$MCMUMU_OK" == true ]] && COMBINED_INPUTS+=("${OUTDIR}/mc_mumugamma.root")

if [[ "${#COMBINED_INPUTS[@]}" -ge 2 ]]; then
    hadd -f "${OUTDIR}/all_merged.root" "${COMBINED_INPUTS[@]}"
    echo "  Size: $(du -h "${OUTDIR}/all_merged.root" | cut -f1)"
elif [[ "${#COMBINED_INPUTS[@]}" -eq 1 ]]; then
    cp "${COMBINED_INPUTS[0]}" "${OUTDIR}/all_merged.root"
    echo -e "${YELLOW}  Warning: only 1 category available, all_merged.root = copy of that.${NC}"
else
    echo -e "${RED}  No category files available — cannot build combined file.${NC}"
fi
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo -e "${YELLOW}============================================================${NC}"
echo -e "${GREEN}  Done${NC}"
echo -e "${YELLOW}============================================================${NC}"
echo ""
echo "Output files in ${OUTDIR}/"
for f in egam3.root egam4.root mc_eegamma.root mc_mumugamma.root all_merged.root; do
    fpath="${OUTDIR}/${f}"
    if [[ -f "${fpath}" ]]; then
        printf "  %-25s  %s\n" "$f" "$(du -h "$fpath" | cut -f1)"
    else
        printf "  %-25s  MISSING\n" "$f"
    fi
done
echo ""
