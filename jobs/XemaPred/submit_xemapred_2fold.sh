#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit held-out within-cohort 2-fold informative-prior jobs for
# PatientSpecificXemaPredV2.
#
# Cross-fitting:
#
#   Target Fold 1 <- Population prior from Fold 2
#   Target Fold 2 <- Population prior from Fold 1
#
# Usage:
#
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_2fold_prior.sh PFDC 1
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_2fold_prior.sh PFDC 4
#
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_2fold_prior.sh Derexyl 1
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_2fold_prior.sh Derexyl 4
#
# Per dataset/horizon:
#
#   2 target folds x 9 PO-SCORAD items = 18 jobs
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"
HORIZON="${2:-1}"

PBS_SCRIPT="jobs/PatientSpecificXemaPred/run_patient_specific_2fold_prior.pbs"

CONFIG_ROOT="configs/${DATASET}/PatientSpecificXemaPred-2foldPrior/t_horizon_${HORIZON}"

TARGET_FOLDS=(1 2)

EXTENT_SIGN="extent"
SUBJ_SIGNS=(itching sleep)
INTENSITY_SIGNS=(dryness redness swelling oozing thickening scratching)

FAMILY_EXTENT="binmc"
FAMILY_SUBJ="binrw"
FAMILY_INT="orderedrw"


# -----------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------

case "$DATASET" in
  PFDC|Derexyl)
    ;;
  *)
    echo "[ERROR] DATASET must be PFDC or Derexyl."
    exit 1
    ;;
esac


case "$HORIZON" in
  1|4)
    ;;
  *)
    echo "[ERROR] HORIZON must be 1 or 4."
    exit 1
    ;;
esac


if [[ ! -f "$PBS_SCRIPT" ]]; then
  echo "[ERROR] Missing PBS script:"
  echo "        $PBS_SCRIPT"
  exit 1
fi


if [[ ! -d "$CONFIG_ROOT" ]]; then
  echo "[ERROR] Missing config root:"
  echo "        $CONFIG_ROOT"
  exit 1
fi


# -----------------------------------------------------------------------
# Submit one PatientSpecificXemaPred job
# -----------------------------------------------------------------------

submit_one () {

  local target_fold="$1"
  local sign="$2"
  local family="$3"

  local source_fold

  if [[ "$target_fold" -eq 1 ]]; then
    source_fold=2
  else
    source_fold=1
  fi


  local ds_lower
  ds_lower="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"


  local CONFIG_YAML
  CONFIG_YAML="${CONFIG_ROOT}/fold_${target_fold}/${ds_lower}_${sign}_${family}.yaml"


  if [[ ! -f "$CONFIG_YAML" ]]; then
    echo "[WARN] Missing config:"
    echo "       $CONFIG_YAML"
    return 0
  fi


  # ---------------------------------------------------------------------
  # Safety check: ensure the generated config has the intended folds.
  # ---------------------------------------------------------------------

  local cfg_target
  local cfg_source

  cfg_target="$(
    awk '
      /^[[:space:]]*patient_fold:[[:space:]]*/ {
        gsub(/#.*/, "")
        sub(/^[[:space:]]*patient_fold:[[:space:]]*/, "")
        gsub(/[[:space:]]/, "")
        print
        exit
      }
    ' "$CONFIG_YAML"
  )"

  cfg_source="$(
    awk '
      /^[[:space:]]*prior_source_fold:[[:space:]]*/ {
        gsub(/#.*/, "")
        sub(/^[[:space:]]*prior_source_fold:[[:space:]]*/, "")
        gsub(/[[:space:]]/, "")
        print
        exit
      }
    ' "$CONFIG_YAML"
  )"


  if [[ "$cfg_target" != "$target_fold" ]]; then
    echo "[ERROR] Target-fold mismatch:"
    echo "        Config:   $CONFIG_YAML"
    echo "        Expected: $target_fold"
    echo "        Found:    $cfg_target"
    exit 1
  fi


  if [[ "$cfg_source" != "$source_fold" ]]; then
    echo "[ERROR] Source-fold mismatch:"
    echo "        Config:   $CONFIG_YAML"
    echo "        Expected: $source_fold"
    echo "        Found:    $cfg_source"
    exit 1
  fi


  if [[ "$cfg_target" == "$cfg_source" ]]; then
    echo "[ERROR] Leakage detected:"
    echo "        patient_fold      = $cfg_target"
    echo "        prior_source_fold = $cfg_source"
    echo "        Config: $CONFIG_YAML"
    exit 1
  fi


  # ---------------------------------------------------------------------
  # Logs
  # ---------------------------------------------------------------------

  local LOG_ROOT
  LOG_ROOT="logs/PatientSpecificXemaPred-2foldPrior/${DATASET}/t_horizon_${HORIZON}/target_fold_${target_fold}/${sign}"

  mkdir -p "$LOG_ROOT"


  # ---------------------------------------------------------------------
  # Job name
  # ---------------------------------------------------------------------

  local jobname

  jobname="$(printf "%s_ps_h%s_t%s_s%s_%s" \
    "$ds_lower" \
    "$HORIZON" \
    "$target_fold" \
    "$source_fold" \
    "$sign"
  )"


  # ---------------------------------------------------------------------
  # Submit
  # ---------------------------------------------------------------------

  qsub \
    -N "$jobname" \
    -o "${LOG_ROOT}/job.${jobname}.out" \
    -e "${LOG_ROOT}/job.${jobname}.err" \
    -v CONFIG_YAML="$CONFIG_YAML",DATASET="$DATASET",SIGN="$sign",HORIZON="$HORIZON",TARGET_FOLD="$target_fold",SOURCE_FOLD="$source_fold" \
    "$PBS_SCRIPT"


  echo "[INFO] Submitted $jobname"
  echo "       Target fold: $target_fold"
  echo "       Source fold: $source_fold"
  echo "       Score:       $sign"
  echo "       Config:      $CONFIG_YAML"
  echo "       Logs:"
  echo "         ${LOG_ROOT}/job.${jobname}.out"
  echo "         ${LOG_ROOT}/job.${jobname}.err"

  sleep 2
}


# -----------------------------------------------------------------------
# Submit all jobs
# -----------------------------------------------------------------------

echo "============================================================"
echo "Submitting PatientSpecificXemaPred cross-fold prior jobs"
echo
echo "DATASET:     $DATASET"
echo "HORIZON:     $HORIZON"
echo "CONFIG_ROOT: $CONFIG_ROOT"
echo "PBS_SCRIPT:  $PBS_SCRIPT"
echo
echo "Cross-fitting:"
echo "  Target Fold 1 <- Source Fold 2"
echo "  Target Fold 2 <- Source Fold 1"
echo
echo "Expected jobs: 18"
echo "============================================================"


for target_fold in "${TARGET_FOLDS[@]}"; do

  if [[ "$target_fold" -eq 1 ]]; then
    source_fold=2
  else
    source_fold=1
  fi

  echo
  echo "------------------------------------------------------------"
  echo "Target Fold ${target_fold} <- Population Source Fold ${source_fold}"
  echo "------------------------------------------------------------"


  # Extent
  submit_one \
    "$target_fold" \
    "$EXTENT_SIGN" \
    "$FAMILY_EXTENT"


  # Subjective scores
  for sign in "${SUBJ_SIGNS[@]}"; do

    submit_one \
      "$target_fold" \
      "$sign" \
      "$FAMILY_SUBJ"

  done


  # Intensity signs
  for sign in "${INTENSITY_SIGNS[@]}"; do

    submit_one \
      "$target_fold" \
      "$sign" \
      "$FAMILY_INT"

  done

done


echo
echo "============================================================"
echo "[OK] All available PatientSpecificXemaPred jobs submitted."
echo "============================================================"