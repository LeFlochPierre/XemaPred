#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------
# Submit cross-cohort PatientSpecificXemaPred jobs
#
# Submits BOTH:
#   t_horizon_1
#   t_horizon_4
#
# Expected config names:
#
#   PFDC_dryness_CrossCohort.yaml
#   PFDC_extent_CrossCohort.yaml
#   ...
#
# No weight-specific configs are used.
#
# Usage:
#
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_xemapred_cross_cohort.sh PFDC
#
# or:
#
#   bash jobs/PatientSpecificXemaPred/submit_patient_specific_xemapred_cross_cohort.sh Derexyl
#
# Expected config directories:
#
#   configs/<DATASET>/PatientSpecificXemaPred-cross-cohort/t_horizon_1/
#   configs/<DATASET>/PatientSpecificXemaPred-cross-cohort/t_horizon_4/
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------

DATASET="${1:-PFDC}"

HORIZONS=(1 4)

PBS_SCRIPT="jobs/PatientSpecificXemaPred/run_patient_specific_xemapred.pbs"

BASE_CONFIG_ROOT="configs/${DATASET}/PatientSpecificXemaPred-cross-cohort"

BASE_LOG_ROOT="logs/PatientSpecificXemaPred-cross-cohort/${DATASET}"


# -----------------------------------------------------------------------
# Expected items
# -----------------------------------------------------------------------

ITEMS=(
    "extent"
    "dryness"
    "redness"
    "swelling"
    "oozing"
    "thickening"
    "scratching"
    "itching"
    "sleep"
)


# -----------------------------------------------------------------------
# Validate PBS runner
# -----------------------------------------------------------------------

if [[ ! -f "$PBS_SCRIPT" ]]; then
    echo "[ERROR] PBS script does not exist:"
    echo "  $PBS_SCRIPT"
    exit 1
fi


# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------

echo "============================================================"
echo "Submitting cross-cohort PatientSpecificXemaPred jobs"
echo
echo "DATASET:     $DATASET"
echo "HORIZONS:    ${HORIZONS[*]}"
echo "PBS_SCRIPT:  $PBS_SCRIPT"
echo "CONFIG ROOT: $BASE_CONFIG_ROOT"
echo "============================================================"


# -----------------------------------------------------------------------
# Submit
# -----------------------------------------------------------------------

N_SUBMITTED=0

for HORIZON in "${HORIZONS[@]}"; do

    CONFIG_ROOT="${BASE_CONFIG_ROOT}/t_horizon_${HORIZON}"
    LOG_ROOT="${BASE_LOG_ROOT}/t_horizon_${HORIZON}"

    echo
    echo "------------------------------------------------------------"
    echo "HORIZON: H${HORIZON}"
    echo "CONFIG_ROOT: $CONFIG_ROOT"
    echo "------------------------------------------------------------"


    # -------------------------------------------------------------------
    # Validate horizon config directory
    # -------------------------------------------------------------------

    if [[ ! -d "$CONFIG_ROOT" ]]; then
        echo "[ERROR] Config directory does not exist:"
        echo "  $CONFIG_ROOT"
        exit 1
    fi

    mkdir -p "$LOG_ROOT"


    # -------------------------------------------------------------------
    # Submit exactly one config per item
    # -------------------------------------------------------------------

    for SCORE in "${ITEMS[@]}"; do

        CONFIG_YAML="${CONFIG_ROOT}/${DATASET}_${SCORE}_CrossCohort.yaml"

        if [[ ! -f "$CONFIG_YAML" ]]; then
            echo "[ERROR] Missing config:"
            echo "  $CONFIG_YAML"
            exit 1
        fi


        # ---------------------------------------------------------------
        # PBS job name
        #
        # Examples:
        #   pfdc_cc_h1_extent
        #   pfdc_cc_h4_dryness
        # ---------------------------------------------------------------

        ds_lower="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

        jobname="${ds_lower}_cc_h${HORIZON}_${SCORE}"


        # ---------------------------------------------------------------
        # Logs
        # ---------------------------------------------------------------

        JOB_LOG_DIR="${LOG_ROOT}/${SCORE}"

        mkdir -p "$JOB_LOG_DIR"


        # ---------------------------------------------------------------
        # Submit
        # ---------------------------------------------------------------

        qsub \
            -N "$jobname" \
            -o "${JOB_LOG_DIR}/job.${jobname}.out" \
            -e "${JOB_LOG_DIR}/job.${jobname}.err" \
            -v CONFIG_YAML="$CONFIG_YAML" \
            "$PBS_SCRIPT"


        echo "[INFO] Submitted: $jobname"
        echo "       Config: $CONFIG_YAML"
        echo "       Logs:   $JOB_LOG_DIR"

        N_SUBMITTED=$((N_SUBMITTED + 1))

        sleep 2

    done

done


# -----------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------

echo
echo "============================================================"
echo "Finished submitting cross-cohort jobs."
echo
echo "Dataset:   $DATASET"
echo "Horizons:  ${HORIZONS[*]}"
echo "Submitted: $N_SUBMITTED"
echo "Expected:  $(( ${#HORIZONS[@]} * ${#ITEMS[@]} ))"
echo "============================================================"