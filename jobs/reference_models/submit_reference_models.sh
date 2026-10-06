#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit all reference-model YAML jobs for one dataset
#
# Run:
# bash jobs/reference_models/submit_reference_models.sh PFDC
# bash jobs/reference_models/submit_reference_models.sh Derexyl
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"

PBS_SCRIPT="jobs/reference_models/run_reference_model.pbs"
CONFIG_ROOT="configs/${DATASET}/reference_models"

# Submit H4 configs. Add 1 if you also want H1.
HORIZONS=(1 4)

if [[ ! -d "$CONFIG_ROOT" ]]; then
  echo "[ERROR] Config root does not exist: $CONFIG_ROOT"
  exit 1
fi

if [[ ! -f "$PBS_SCRIPT" ]]; then
  echo "[ERROR] PBS script does not exist: $PBS_SCRIPT"
  exit 1
fi

submit_config () {
  local CONFIG_YAML="$1"

  local ds_lower
  ds_lower="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

  local rel
  rel="${CONFIG_YAML#${CONFIG_ROOT}/}"

  local model
  model="$(echo "$rel" | cut -d "/" -f 1)"

  local horizon_dir
  horizon_dir="$(echo "$rel" | cut -d "/" -f 2)"

  local file
  file="$(basename "$CONFIG_YAML" .yaml)"

  local LOG_ROOT="logs/reference_models/${DATASET}/${model}/${horizon_dir}/${file}"
  mkdir -p "$LOG_ROOT"

  local jobname
  jobname="$(printf "%s_ref_%s_%s" "$ds_lower" "$model" "$file")"

  # Some PBS systems have short job-name limits.
  jobname="${jobname:0:60}"

  qsub -N "$jobname" \
    -o "${LOG_ROOT}/job.${jobname}.out" \
    -e "${LOG_ROOT}/job.${jobname}.err" \
    -v CONFIG_YAML="$CONFIG_YAML" \
    "$PBS_SCRIPT"

  echo "[INFO] Submitted $jobname"
  echo "  CONFIG_YAML: $CONFIG_YAML"
  echo "  Logs: ${LOG_ROOT}/job.${jobname}.out"
}

for h in "${HORIZONS[@]}"; do
  horizon_dir="t_horizon_${h}"

  find "$CONFIG_ROOT" -path "*/${horizon_dir}/*.yaml" -type f | sort | while read -r cfg; do
    submit_config "$cfg"
  done
done