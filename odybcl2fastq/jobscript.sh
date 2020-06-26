#!/bin/bash
# properties = {properties}

singularity exec "$SINGULARITY_CONTAINER" bash -c '{exec_job}'
