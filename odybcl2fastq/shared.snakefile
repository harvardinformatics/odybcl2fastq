'''
snakemake shared rules

Created on  2020-04-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@author: Nathan Weeks <nweeks@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
from odybcl2fastq.bauer_db import BauerDB
from odybcl2fastq.status_db import StatusDB
import odybcl2fastq.util as util
import logging
import os
import shutil
from pathlib import Path

os.umask(0o002) # created files/directories default to mode 0775

# allow an empty suffix
wildcard_constraints:
    suffix=".*"

analysis_dir = Path('/sequencing', 'analysis', config['run'] + config['suffix'])
# set the output_run and status_dir which may include a suffix or a mask_suffix
status_dir_root = 'status_test' if ody_config.TEST else 'status'
status_dir = status_dir_root # default status_dir
if 'mask_suffix' in config and config['mask_suffix']:
    sample_sheet_path = "/sequencing/source/%s/SampleSheet.csv" % (config['run'])
    # for a mask suffix we need a subdir for status for each mask
    status_dir += '/%s' % config['mask_suffix']
else:
    sample_sheet_path = "/sequencing/source/%s/SampleSheet%s.csv" % (config['run'], config['suffix'])

# set up bauer db for step updates
bauer = BauerDB(sample_sheet_path)

onstart:
    """
    touch processed file to prevent reprocessing
    prepare log, script and status dirs
    """
    shell("mkdir -p /sequencing/source/{run}/{status}", run=config['run'], status=status_dir)
    shell("touch /sequencing/source/{run}/{status_root}/ody.processed", run=config['run'], status_root=status_dir_root)
    shell("touch /sequencing/source/{run}/{status}/ody.processed", run=config['run'], status=status_dir)
    update_analysis({'status': 'processing'})

rule insert_run_into_bauer_db:
    """
    insert the run into the bauer_db
    """
    input:
        sample_sheet_path
    output:
        f"/sequencing/source/{config['run']}/{status_dir}/analysis_id"
    run:
        if ody_config.TEST:
            analysis_id = 'test'
        else:
            bauer = BauerDB(sample_sheet_path)
            bauer.insert_run()
            # TODO: consider the implications of storing output_run in the db
            analysis_id = bauer.send_data('requests', {"run": config['run'], "status":"processing", "step":"demultiplex"})
        analysis_file_path = '/sequencing/source/%s/%s/analysis_id' % (config['run'], status_dir)
        with open(analysis_file_path, 'w+') as f:
            f.write(str(analysis_id))
        # for a new analysis remove any existing analysis_dir to ensure a total restart
        if analysis_dir.exists(): shutil.rmtree(analysis_dir)
        analysis_dir.mkdir()
        Path(analysis_dir, 'log').mkdir()
        Path(analysis_dir, 'script').mkdir()

rule update_lims_db:
    """
    connect the submission with the run in lims
    """
    input:
        expand("/sequencing/source/{{run}}/{status}/demultiplex.processed", status=status_dir),
        sample_sheet_path
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/update_lims_db.processed", status=status_dir))
    run:
        if not ody_config.TEST:
            sample_sheet = SampleSheet(sample_sheet_path)
            instrument = sample_sheet.get_instrument()
            update_lims_db(config['run'], sample_sheet.sections, instrument)
            # also update step to fastqc
            update_analysis({'step': 'quality', 'status': 'processing'})

rule fastqc_cmd:
    """
    build a bash file with the fastqc cmd
    """
    input:
        ancient(expand("/sequencing/source/{run}/{status}/demultiplex.processed", run=config['run'], status=status_dir))
    output:
        expand("/sequencing/analysis/{{run}}{{suffix}}/script/fastqc.sh")
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="mkdir -p /sequencing/analysis/{config[run]}{config[suffix]}/QC\n"
        cmd+="cd /sequencing/analysis/{config[run]}{config[suffix]}/fastq/\n"
        cmd+="find . -name '*.fastq.gz' ! -name 'Undetermined*' -exec /usr/bin/time -v fastqc -o /sequencing/analysis/{config[run]}{config[suffix]}/QC --threads \$SLURM_JOB_CPUS_PER_NODE {{}} +"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule fastqc:
    """
    run bash file for fastqc
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        expand("/sequencing/analysis/{{run}}{suffix}/script/fastqc.sh", suffix=config['suffix'])
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/fastqc.processed", status=status_dir))
    shell:
        """
        {input}
        """

rule multiqc:
    """
    run multiqc
    """
    input:
        ancient(expand("/sequencing/source/{run}/{status}/fastqc.processed", run=config['run'], status=status_dir))
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/multiqc.processed", status=status_dir))
    shell:
        """
        cd /sequencing/analysis/{config[run]}{config[suffix]}/QC
        /usr/bin/time -v multiqc /sequencing/analysis/{config[run]}{config[suffix]}/QC
        """

rule cp_source_to_output:
    """
    copy a few files from source to output dir
    """
    input:
        expand("/sequencing/source/{run}/{status}/demultiplex.processed", run=config['run'], status=status_dir)
    params:
        sample_sheet="SampleSheet%s.csv" % config['suffix'],
        run_info="RunInfo.xml",
        interop="InterOp",
        nextseq_run_params="RunParameters.xml",
        hiseq_run_params="runParameters.xml"
    output:
        sample_sheet=expand("/sequencing/analysis/{{run}}{{suffix}}/SampleSheet.csv"),
        run_info=expand("/sequencing/analysis/{{run}}{{suffix}}/RunInfo.xml")
    shell:
        """
        cp /sequencing/source/{config[run]}/{params.sample_sheet} {output.sample_sheet}
        cp /sequencing/source/{config[run]}/{params.run_info} {output.run_info}
        rsync --info=STATS -rtl --safe-links --perms --chmod=Dug=rwx,Fug=rw /sequencing/source/{config[run]}/{params.interop}/ /sequencing/analysis/{config[run]}{config[suffix]}/InterOp/
        # copy these if they exist
        cp /sequencing/source/{config[run]}/{params.nextseq_run_params} /sequencing/analysis/{config[run]}{config[suffix]}/{params.nextseq_run_params} 2>/dev/null || :
        cp /sequencing/source/{config[run]}/{params.hiseq_run_params} /sequencing/analysis/{config[run]}{config[suffix]}/{params.hiseq_run_params} 2>/dev/null || :

        """

rule checksum:
    """
    calculate checksum for all the fastq files
    """
    input:
        expand("/sequencing/source/{run}/{status}/demultiplex.processed", run=config['run'], status=status_dir)
    output:
        checksum=expand("/sequencing/analysis/{{run}}{{suffix}}/md5sum.txt"),
    shell:
        """
        files=$(find /sequencing/analysis/{config[run]}{config[suffix]}/ -name *.fastq.gz -print0 | xargs -0)
        md5sum $files > {output.checksum}
        """

def update_analysis(data):
    if not ody_config.TEST:
        analysis_file_path = '/sequencing/source/%s/%s/analysis_id' % (config['run'], status_dir)
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
