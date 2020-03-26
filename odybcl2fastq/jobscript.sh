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

# check if singularity img exist
if ! [ -f $ODY_SING_IMG ]; then
    echo "file does not exist: $ODY_SING_IMG"
    exit 1
    break
fi

singularity exec -B $ODY_SOURCE:/source \
                -B $published:/published \
                -B $analysis:/output \
                -B /n/home_rc/mportermahoney/repos/odybcl2fastq_10x_dev/odybcl2fastq:/app \
                -B $ODY_SNAKEMAKE_WORKDIR:/snakemake \
                $ODY_SING_IMG bash -c '{exec_job}'
