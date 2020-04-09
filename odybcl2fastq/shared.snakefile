'''
snakemake shared rules

Created on  2020-04-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@author: Nathan Weeks <nweeks@fas.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
from odybcl2fastq.bauer_db import BauerDB
import odybcl2fastq.util as util
import pandas as pd
import os

# set the output_run and status_dir which may include a suffix or a mask_suffix
status_dir_root = 'status_test' if ody_config.TEST else 'status'
status_dir = status_dir_root # default status_dir
output_run = config['output_run']
sample_sheet_name = config['sample_sheet_name']
sample_sheet_path = "%s%s/%s" % (ody_config.SOURCE_DIR, config['run'], sample_sheet_name)

# set up bauer db for step updates
bauer = BauerDB(sample_sheet_path)

# allow an empty suffix
wildcard_constraints:
    run_suffix=".*"

onstart:
    """
    touch processed file to prevent reprocessing
    prepare log, script and status dirs
    """
    shell("mkdir -p {source}{run}/{status}", source={ody_config.SOURCE_DIR}, run=config['run'], status=status_dir)
    shell("touch {source}{run}/{status_root}/ody.processed", source=ody_config.SOURCE_DIR, run=config['run'], status_root=status_dir_root)
    shell("touch {source}{run}/{status}/ody.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)
    shell("mkdir -p {output}{run}", output={ody_config.OUTPUT_DIR}, run=output_run)
    shell("mkdir -p {output}{run}/log", output={ody_config.OUTPUT_DIR}, run=output_run)
    shell("mkdir -p {output}{run}/script", output={ody_config.OUTPUT_DIR}, run=output_run)
    update_analysis({'status': 'processing'})

rule insert_run_into_bauer_db:
    """
    insert the run into the bauer_db
    """
    input:
        sample_sheet=expand("{source}{{run}}/{sample_sheet_name}", source=ody_config.SOURCE_DIR, sample_sheet_name=sample_sheet_name)
    output:
        expand("{source}{{run}}/{status}/analysis_id", source=ody_config.SOURCE_DIR, status=status_dir)
    run:
        if ody_config.TEST:
            analysis_id = 'test'
        else:
            bauer = BauerDB(input.sample_sheet[0])
            bauer.insert_run()
            # TODO: consider the implications of storing output_run in the db
            analysis_id = bauer.send_data('requests', {"run": output_run, "status":"processing", "step":"demultiplex"})
        analysis_file_path = '%s%s/%s/analysis_id' % (ody_config.SOURCE_DIR, wildcards.run, status_dir)
        with open(analysis_file_path, 'w+') as f:
            f.write(str(analysis_id))
        # for a new analysis remove the script dir to ensure a total restart
        shell("rm -f {output}{run}/script/*", output={ody_config.OUTPUT_DIR}, run=output_run)

rule update_lims_db:
    """
    connect the submission with the run in lims
    """
    input:
        expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=status_dir),
        sample_sheet=expand("{source}{{run}}/{sample_sheet_name}", source=ody_config.SOURCE_DIR, sample_sheet_name=sample_sheet_name)
    output:
        touch(expand("{source}{{run}}/{status}/update_lims_db.processed", source=ody_config.SOURCE_DIR, status=status_dir))
    run:
        if not ody_config.TEST:
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
        ancient(expand("{source}{run}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir))
    output:
        expand("{output}{{output_run}}/script/fastqc.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/QC\n"
        cmd+="cd {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/fastq/\n"
        cmd+="find . -name '*.fastq.gz' ! -name 'Undetermined*' -exec /usr/bin/time -v fastqc -o {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/QC --threads \$SLURM_JOB_CPUS_PER_NODE {{}} +"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule fastqc:
    """
    run bash file for fastqc
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        expand("{output}{output_run}/script/fastqc.sh", output=ody_config.OUTPUT_DIR, output_run=output_run)
    output:
        touch(expand("{source}{{run}}/{status}/fastqc.processed", source=ody_config.SOURCE_DIR, status=status_dir))
    shell:
        """
        {input}
        """

rule cp_source_to_output:
    """
    copy a few files from source to output dir
    """
    input:
        expand("{source}{run}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)
    params:
        sample_sheet=sample_sheet_name,
        run_info="RunInfo.xml",
        interop="InterOp",
        nextseq_run_params="RunParameters.xml",
        hiseq_run_params="runParameters.xml"
    output:
        sample_sheet=expand("{output}{{output_run}}/SampleSheet.csv", output=ody_config.OUTPUT_DIR),
        run_info=expand("{output}{{output_run}}/RunInfo.xml", output=ody_config.OUTPUT_DIR)
    shell:
        """
        cp {ody_config.SOURCE_DIR}{config[run]}/{params.sample_sheet} {output.sample_sheet}
        cp {ody_config.SOURCE_DIR}{config[run]}/{params.run_info} {output.run_info}
        rsync --info=STATS -rtl --safe-links --perms --chmod=Dug=rwx,Fug=rw {ody_config.SOURCE_DIR}{config[run]}/{params.interop}/ {ody_config.OUTPUT_DIR}{wildcards.output_run}/InterOp/
        # copy these if they exist
        cp {ody_config.SOURCE_DIR}{config[run]}/{params.nextseq_run_params} {ody_config.OUTPUT_DIR}{wildcards.output_run}/{params.nextseq_run_params} 2>/dev/null || :
        cp {ody_config.SOURCE_DIR}{config[run]}/{params.hiseq_run_params} {ody_config.OUTPUT_DIR}{wildcards.output_run}/{params.hiseq_run_params} 2>/dev/null || :

        """

rule checksum:
    """
    calculate checksum for all the fastq files
    """
    input:
        expand("{source}{run}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)
    output:
        checksum=expand("{output}{{output_run}}/md5sum.txt", output=ody_config.OUTPUT_DIR),
    shell:
        """
        files=$(find {ody_config.OUTPUT_DIR}{wildcards.output_run}/ -name *.fastq.gz -print0 | xargs -0)
        md5sum $files > {output.checksum}
        """

def update_analysis(data):
    if not ody_config.TEST:
        analysis_file_path = '%s%s/%s/analysis_id' % (ody_config.SOURCE_DIR, config['run'], status_dir)
        if os.path.isfile(analysis_file_path):
            with open(analysis_file_path, 'r') as ln:
                analysis_id = ln.readline().strip()
            bauer.update_data('requests', analysis_id, data)

def get_submissions(sample_sheet, instrument):
    subs = set()
    for key, row in sample_sheet['Data'].items():
        if row['Description']:
            subs.add(row['Description'])
    return list(subs)

def update_lims_db(run, sample_sheet, instrument):
    runlogger = logging.getLogger('run_logger')
    runlogger.info('Start db update for %s\n' % run)
    subs = get_submissions(sample_sheet, instrument)
    stdb = StatusDB()
    stdb.link_run_and_subs(run, subs)
    analysis = stdb.insert_analysis(run, ', '.join(subs))
    runlogger.info('End db update for %s\n' % analysis)
