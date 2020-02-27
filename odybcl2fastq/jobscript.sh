#!/bin/bash
# properties = {properties}
singularity exec -B $ODY_SOURCE:/source \
                -B $ODY_SEQ_ROOT/published:/published \
                -B $ODY_SEQ_ROOT/analysis:/output \
                -B $ODY_SNAKEMAKE_WORKDIR:/snakemake \
                $ODY_SING_IMG bash -c '{exec_job}'
