'''
snakemake shared rules

Created on  2020-04-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

from odybcl2fastq.run import update_lims_db, setup_run_logger
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
from odybcl2fastq.bauer_db import BauerDB
import odybcl2fastq.util as util
import pandas as pd
import os

STATUS_DIR = 'status_test' if ody_config.TEST else 'status'

# set up bauer db for step updates
sample_sheet_path = "%s%s/SampleSheet.csv" % (ody_config.SOURCE_DIR, config['run'])
bauer = BauerDB(sample_sheet_path)

rule insert_run_into_bauer_db:
    """
    insert the run into the bauer_db
    """
    input:
        sample_sheet=expand("{source}{{run}}/SampleSheet.csv", source=ody_config.SOURCE_DIR)
    output:
        expand("{source}{{run}}/{status}/analysis_id", source=ody_config.SOURCE_DIR, status=STATUS_DIR)
    run:
        if ody_config.TEST:
            analysis_id = 'test'
        else:
            bauer = BauerDB(input.sample_sheet[0])
            bauer.insert_run()
            analysis_id = bauer.send_data('requests', {"run": wildcards.run, "status":"processing", "step":"demultiplex"})
        analysis_file_path = '%s%s/%s/analysis_id' % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR)
        with open(analysis_file_path, 'w+') as f:
            f.write(str(analysis_id))
        # for a new analysis remove the script dir to ensure a total restart
        shell("rm -f {output}{run}/script/*", output={ody_config.OUTPUT_DIR}, run=config['run'])

rule update_lims_db:
    """
    connect the submission with the run in lims
    """
    input:
        expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        sample_sheet=expand("{source}{{run}}/SampleSheet.csv", source=ody_config.SOURCE_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/update_lims_db.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    run:
        if not ody_config.TEST:
            setup_run_logger(wildcards.run, False)
            sample_sheet = SampleSheet(input.sample_sheet[0])
            instrument = sample_sheet.get_instrument()
            run_folder = ody_config.SOURCE_DIR + wildcards.run
            update_lims_db(run_folder, sample_sheet.sections, instrument)
            # also update step to fastqc
            update_analysis({'step': 'quality', 'status': 'processing'})

rule fastqc_cmd:
    """
    build a bash file with the fastqc cmd
    """
    input:
        ancient(expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    output:
        expand("{output}{{run}}/script/fastqc.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/QC\n"
        cmd+="cd {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq/\n"
        cmd+="find . -name '*.fastq.gz' ! -name 'Undetermined*' -exec /usr/bin/time -v fastqc -o {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/QC --threads \$SLURM_JOB_CPUS_PER_NODE {{}} +"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule fastqc:
    """
    run bash file for fastqc
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        expand("{output}{{run}}/script/fastqc.sh", output=ody_config.OUTPUT_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/fastqc.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    shell:
        """
        {input}
        """

rule cp_source_to_output:
    """
    copy a few files from source to output dir
    """
    input:
        expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR)
    params:
        sample_sheet="SampleSheet.csv",
        run_info="RunInfo.xml",
        interop="InterOp",
        nextseq_run_params="RunParameters.xml",
        hiseq_run_params="runParameters.xml"
    output:
        sample_sheet=expand("{output}{{run}}/SampleSheet.csv", output=ody_config.OUTPUT_DIR),
        run_info=expand("{output}{{run}}/RunInfo.xml", output=ody_config.OUTPUT_DIR)
    shell:
        """
        cp {ody_config.SOURCE_DIR}{wildcards.run}/{params.sample_sheet} {output.sample_sheet}
        cp {ody_config.SOURCE_DIR}{wildcards.run}/{params.run_info} {output.run_info}
        rsync --info=STATS -rtl --safe-links --perms --chmod=Dug=rwx,Fug=rw {ody_config.SOURCE_DIR}{wildcards.run}/{params.interop}/ {ody_config.OUTPUT_DIR}{wildcards.run}/InterOp/
        # copy these if they exist
        cp {ody_config.SOURCE_DIR}{wildcards.run}/{params.nextseq_run_params} {ody_config.OUTPUT_DIR}{wildcards.run}/{params.nextseq_run_params} 2>/dev/null || :
        cp {ody_config.SOURCE_DIR}{wildcards.run}/{params.hiseq_run_params} {ody_config.OUTPUT_DIR}{wildcards.run}/{params.hiseq_run_params} 2>/dev/null || :

        """

rule checksum:
    """
    calculate checksum for all the fastq files
    """
    input:
        expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR)
    output:
        checksum=expand("{output}{{run}}/md5sum.txt", output=ody_config.OUTPUT_DIR),
    shell:
        """
        files=$(find {ody_config.OUTPUT_DIR}{wildcards.run}/ -name *.fastq.gz -print0 | xargs -0)
        md5sum $files > {output.checksum}
        """

def update_analysis(data):
    if not ody_config.TEST:
        analysis_file_path = '%s%s/%s/analysis_id' % (ody_config.SOURCE_DIR, config['run'], STATUS_DIR)
        if os.path.isfile(analysis_file_path):
            with open(analysis_file_path, 'r') as ln:
                analysis_id = ln.readline().strip()
            bauer.update_data('requests', analysis_id, data)
