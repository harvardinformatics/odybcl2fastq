#!/bin/bash
# properties = {properties}
published="$ODY_SEQ_ROOT/published"
analysis="$ODY_SEQ_ROOT/analysis"

# check if singularity mount paths exist
for path in \
   $ODY_SOURCE \
   $analysis \
   $published \
   $ODY_SNAKEMAKE_WORKDIR
do
    if ! [ -d $path ]; then
        echo "path does not exist: $path"
        exit 1
        break
    fi
done

singularity exec -B $ODY_SOURCE:/source \
                -B $published:/published \
                -B $analysis:/output \
                -B $ODY_SNAKEMAKE_WORKDIR:/snakemake \
                "${SINGULARITY_CONTAINER}" bash -c '{exec_job}'
