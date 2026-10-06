#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit all EczemaPred item-model jobs
#
# Models:
# - extent: BinMC
# - intensity signs: OrderedRW
# - subjective symptoms: BinRW
#
# Run:
#   bash jobs/EczemaPred/submit_eczemapred.sh PFDC
#   bash jobs/EczemaPred/submit_eczemapred.sh Derexyl
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"

PBS_SCRIPT="jobs/EczemaPred/run_eczemapred.pbs"
CONFIG_ROOT="configs/${DATASET}/EczemaPred"

HORIZONS=(1 4)

EXTENT_SIGN="extent"
SUBJ_SIGNS=(itching sleep)
INTENSITY_SIGNS=(dryness redness swelling oozing thickening scratching)

submit_one () {
  local horizon="$1"
  local sign="$2"
  local family="$3"

  local ds_lower
  ds_lower="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

  local CONFIG_YAML="${CONFIG_ROOT}/t_horizon_${horizon}/${ds_lower}_${sign}_${family}.yaml"

  if [[ ! -f "$CONFIG_YAML" ]]; then
    echo "[WARN] Missing config: $CONFIG_YAML"
    return 0
  fi

  local LOG_ROOT="logs/EczemaPred/${DATASET}/t_horizon_${horizon}/${sign}"
  mkdir -p "$LOG_ROOT"

  local jobname
  jobname="$(printf "%s_eczemapred_h%s_%s" "$ds_lower" "$horizon" "$sign")"

  qsub -N "$jobname" \
    -o "${LOG_ROOT}/job.${jobname}.out" \
    -e "${LOG_ROOT}/job.${jobname}.err" \
    -v CONFIG_YAML="$CONFIG_YAML" \
    "$PBS_SCRIPT"

  echo "[INFO] Submitted $jobname"
  echo "  CONFIG_YAML: $CONFIG_YAML"
  echo "  Logs: ${LOG_ROOT}/job.${jobname}.out"
  echo ""

  sleep 2
}

for h in "${HORIZONS[@]}"; do
  submit_one "$h" "$EXTENT_SIGN" "binmc"

  for s in "${SUBJ_SIGNS[@]}"; do
    submit_one "$h" "$s" "binrw"
  done

  for s in "${INTENSITY_SIGNS[@]}"; do
    submit_one "$h" "$s" "orderedrw"
  done
done