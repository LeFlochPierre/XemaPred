#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit all XemaPred item-model jobs
#
# Usage:
# bash jobs/XemaPred/submit_xemapred.sh PFDC
# bash jobs/XemaPred/submit_xemapred.sh Derexyl
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"
PBS_SCRIPT="jobs/XemaPred/run_xemapred.pbs"

CONFIG_ROOT="configs/${DATASET}/XemaPred"

HORIZONS=(1 4)

EXTENT_SIGN="extent"
SUBJ_SIGNS=(itching sleep)
INTENSITY_SIGNS=(dryness redness swelling oozing thickening scratching)

FAMILY_EXTENT="binmc"
FAMILY_SUBJ="binrw"
FAMILY_INT="orderedrw"

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

  local LOG_ROOT="logs/XemaPred/${DATASET}/t_horizon_${horizon}/${sign}"
  mkdir -p "$LOG_ROOT"

  local jobname
  jobname="$(printf "%s_xemapred_h%s_%s" "$ds_lower" "$horizon" "$sign")"

  qsub -N "$jobname" \
    -o "${LOG_ROOT}/job.${jobname}.out" \
    -e "${LOG_ROOT}/job.${jobname}.err" \
    -v CONFIG_YAML="$CONFIG_YAML",DATASET="$DATASET",SIGN="$sign",HORIZON="$horizon" \
    "$PBS_SCRIPT"

  echo "[INFO] Submitted $jobname"
  echo "  CONFIG_YAML: $CONFIG_YAML"
  echo "  Logs:"
  echo "    ${LOG_ROOT}/job.${jobname}.out"
  echo "    ${LOG_ROOT}/job.${jobname}.err"

  sleep 2
}

echo "======================================"
echo "Submitting XemaPred jobs"
echo "DATASET: $DATASET"
echo "CONFIG_ROOT: $CONFIG_ROOT"
echo "PBS_SCRIPT: $PBS_SCRIPT"
echo "======================================"

for h in "${HORIZONS[@]}"; do
  submit_one "$h" "$EXTENT_SIGN" "$FAMILY_EXTENT"

  for s in "${SUBJ_SIGNS[@]}"; do
    submit_one "$h" "$s" "$FAMILY_SUBJ"
  done

  for s in "${INTENSITY_SIGNS[@]}"; do
    submit_one "$h" "$s" "$FAMILY_INT"
  done
done

echo "[INFO] All available XemaPred jobs submitted."
