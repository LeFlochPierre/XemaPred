#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit one complete 2-fold mixture-weight experiment.
#
# Usage:
#
#   bash jobs/PatientSpecificXemaPred/submit_2fold_mixture.sh PFDC W100
#   bash jobs/PatientSpecificXemaPred/submit_2fold_mixture.sh Derexyl W100
#
# Examples:
#   W000 = 0% population / 100% base
#   W005 = 5% population / 95% base
#   W010 = 10% population / 90% base
#   W020 = 20% population / 80% base
#   W050 = 50% population / 50% base
#   W100 = 100% population
#
# One invocation submits:
#   2 horizons x 2 target folds x 9 PO-SCORAD items = 36 PBS jobs.
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"
WEIGHT_TAG="${2:-W100}"

PBS_SCRIPT="jobs/PatientSpecificXemaPred/run_patient_specific_xemapred_2folds.pbs"
CONFIG_ROOT="configs/${DATASET}/PatientSpecificXemaPred-2fold"

HORIZONS=(1 4)
TARGET_FOLDS=(1 2)

EXTENT_SIGN="extent"
SUBJ_SIGNS=(itching sleep)
INTENSITY_SIGNS=(dryness redness swelling oozing thickening scratching)

FAMILY_EXTENT="binmc"
FAMILY_SUBJ="binrw"
FAMILY_INT="orderedrw"

# -----------------------------------------------------------------------
# Basic validation
# -----------------------------------------------------------------------

if [[ "$DATASET" != "PFDC" && "$DATASET" != "Derexyl" ]]; then
  echo "[ERROR] DATASET must be PFDC or Derexyl."
  exit 1
fi

if [[ ! "$WEIGHT_TAG" =~ ^W[0-9]{3}$ ]]; then
  echo "[ERROR] WEIGHT_TAG must look like W000, W005, W010, W020, W050, W100."
  exit 1
fi

if [[ ! -f "$PBS_SCRIPT" ]]; then
  echo "[ERROR] PBS script does not exist:"
  echo "        $PBS_SCRIPT"
  exit 1
fi

if [[ ! -d "$CONFIG_ROOT" ]]; then
  echo "[ERROR] Config root does not exist:"
  echo "        $CONFIG_ROOT"
  exit 1
fi

ds_lower="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

if [[ "$DATASET" == "PFDC" ]]; then
  ds_short="pf"
else
  ds_short="de"
fi

wcode="${WEIGHT_TAG#W}"

# -----------------------------------------------------------------------
# Submit one item/fold/horizon
# -----------------------------------------------------------------------

submit_one () {
  local horizon="$1"
  local target_fold="$2"
  local source_fold="$3"
  local sign="$4"
  local family="$5"

  local CONFIG_YAML="${CONFIG_ROOT}/t_horizon_${horizon}/fold_${target_fold}/${ds_lower}_${sign}_${family}.yaml"

  if [[ ! -f "$CONFIG_YAML" ]]; then
    echo "[WARN] Missing config: $CONFIG_YAML"
    return 0
  fi

  # Sanity-check the cross-fit direction encoded in the YAML.
  local cfg_target
  local cfg_source
  local cfg_combination

  cfg_target="$(awk '$1=="patient_fold:" {print $2; exit}' "$CONFIG_YAML")"
  cfg_source="$(awk '$1=="prior_source_fold:" {print $2; exit}' "$CONFIG_YAML")"
  cfg_combination="$(awk '$1=="prior_combination:" {print $2; exit}' "$CONFIG_YAML" | tr -d '"'"'")"

  if [[ "$cfg_target" != "$target_fold" ]]; then
    echo "[ERROR] patient_fold mismatch in $CONFIG_YAML"
    echo "        expected=$target_fold found=$cfg_target"
    exit 1
  fi

  if [[ "$cfg_source" != "$source_fold" ]]; then
    echo "[ERROR] prior_source_fold mismatch in $CONFIG_YAML"
    echo "        expected=$source_fold found=$cfg_source"
    exit 1
  fi

  if [[ "$cfg_target" == "$cfg_source" ]]; then
    echo "[ERROR] Leakage detected in $CONFIG_YAML"
    exit 1
  fi

  if [[ "$cfg_combination" != "mixture" ]]; then
    echo "[ERROR] Expected prior_combination: mixture in $CONFIG_YAML"
    echo "        found=$cfg_combination"
    exit 1
  fi

  local LOG_ROOT="logs/PatientSpecificXemaPredMixture/${DATASET}/${WEIGHT_TAG}/t_horizon_${horizon}/fold_${target_fold}/${sign}"
  mkdir -p "$LOG_ROOT"

  # Keep PBS job name compact.
  local sign_short="${sign:0:3}"
  local jobname="${ds_short}${wcode}h${horizon}f${target_fold}${sign_short}"

  qsub -N "$jobname" \
    -o "${LOG_ROOT}/job.${jobname}.out" \
    -e "${LOG_ROOT}/job.${jobname}.err" \
    -v CONFIG_YAML="$CONFIG_YAML",DATASET="$DATASET",SIGN="$sign",TARGET_FOLD="$target_fold",SOURCE_FOLD="$source_fold",HORIZON="$horizon",WEIGHT_TAG="$WEIGHT_TAG" \
    "$PBS_SCRIPT"

  echo "[INFO] Submitted $jobname"
  echo "  Dataset:     $DATASET"
  echo "  Weight:      $WEIGHT_TAG"
  echo "  Horizon:     H$horizon"
  echo "  Target fold: $target_fold"
  echo "  Source fold: $source_fold"
  echo "  Item:        $sign"
  echo "  CONFIG_YAML: $CONFIG_YAML"
  echo "  Logs:        $LOG_ROOT"
  echo

  sleep 1
}

echo "============================================================"
echo "Submitting PatientSpecificXemaPred 2-fold mixture experiment"
echo "DATASET:     $DATASET"
echo "WEIGHT_TAG:  $WEIGHT_TAG"
echo "CONFIG_ROOT: $CONFIG_ROOT"
echo "PBS_SCRIPT:  $PBS_SCRIPT"
echo "============================================================"

n_submitted=0

for h in "${HORIZONS[@]}"; do

  for target_fold in "${TARGET_FOLDS[@]}"; do

    if [[ "$target_fold" -eq 1 ]]; then
      source_fold=2
    else
      source_fold=1
    fi

    submit_one "$h" "$target_fold" "$source_fold" "$EXTENT_SIGN" "$FAMILY_EXTENT"
    n_submitted=$((n_submitted + 1))

    for s in "${SUBJ_SIGNS[@]}"; do
      submit_one "$h" "$target_fold" "$source_fold" "$s" "$FAMILY_SUBJ"
      n_submitted=$((n_submitted + 1))
    done

    for s in "${INTENSITY_SIGNS[@]}"; do
      submit_one "$h" "$target_fold" "$source_fold" "$s" "$FAMILY_INT"
      n_submitted=$((n_submitted + 1))
    done

  done
done

echo "============================================================"
echo "[OK] Submission pass complete."
echo "Requested jobs: $n_submitted"
echo "Dataset:        $DATASET"
echo "Weight:         $WEIGHT_TAG"
echo "============================================================"
