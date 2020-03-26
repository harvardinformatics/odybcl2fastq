'''
snakemake workflow for non 10x

Created on  2020-03-26

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

configfile: "/app/odybcl2fastq/snakemake_non_10x_config.json"
localrules: all, update_lims_db, cp_source_to_output, checksum, publish, demultiplex_cmd, fastqc_cmd, insert_run_into_bauer_db
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq.parsers import parse_stats
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
import odybcl2fastq.util as util
import os

include: "shared.snakefile"

# get some information about the run
sample_sheet = SampleSheet(sample_sheet_path)
instrument = sample_sheet.get_instrument()
run_info = ody_config.SOURCE_DIR + config['run'] + '/RunInfo.xml'
run_type = sample_sheet.get_run_type()

# get basemasks for entire run and then use the one for the indexing strategy
# passed into config['suffix'], for runs with only one strategy suffix will be
# an empty string
mask_lists, mask_samples = extract_basemasks(sample_sheet.sections['Data'], run_info, instrument, run_type)
jobs_tot = len(mask_lists)
if jobs_tot > 1:
    mask = config['suffix'].replace('_', ',')
    sample_sheet_path = sample_sheet.write_new_sample_sheet(mask_samples[mask], config['suffix'])
else:
    mask = next(iter(mask_lists))
mask_opt_list = mask_lists[mask]
mask_switch = '--use-bases-mask'
mask_opt = (mask_switch + ' ' + (' ' + mask_switch + ' ').join(mask_opt_list))

onstart:
    """
    touch processed file to prevent reprocessing
    prepare log, script and status dirs
    """
    shell("mkdir -p {source}{run}/{status}", source={ody_config.SOURCE_DIR}, run=config['run'], status=STATUS_DIR)
    shell("touch {source}{run}/{status}/ody.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=STATUS_DIR)
    shell("mkdir -p {output}{run}", output={ody_config.OUTPUT_DIR}, run=config['run'])
    shell("mkdir -p {output}{run}/log", output={ody_config.OUTPUT_DIR}, run=config['run'])
    shell("mkdir -p {output}{run}/script", output={ody_config.OUTPUT_DIR}, run=config['run'])
    update_analysis({'status': 'processing'})

rule all:
    """
    final output of workflow
    """
    input:
        expand("{source}{run}/{status}/ody.complete", source=ody_config.SOURCE_DIR, run=config['run'], status=STATUS_DIR)

rule demultiplex_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("{source}{{run}}/{status}/analysis_id", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        run_dir=expand("{source}{{run}}/SampleSheet.csv", source=ody_config.SOURCE_DIR)
    output:
        expand("{output}{{run}}/script/demultiplex.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq\n"
        cmd+="/usr/bin/time -v bcl2fastq --adapter-stringency 0.9 --barcode-mismatches 0 --fastq-compression-level 4 --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --min-log-level INFO --minimum-trimmed-read-length 1 --sample-sheet {ody_config.SOURCE_CLUSTER_PATH}{wildcards.run}/SampleSheet.csv --runfolder-dir {ody_config.SOURCE_CLUSTER_PATH}{wildcards.run} --output-dir {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq --processing-threads 8 {mask_opt} || exit_code=\$?\n"
        cmd+="exit \$exit_code"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule demultiplex:
    """
    run bash file for demux
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        expand("{output}{{run}}/script/demultiplex.sh", output=ody_config.OUTPUT_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

def publish_input(wildcards):
    """
    determine which files need to be ready to publish the run
    count is not run if there is not reference genome
    """
    input = {
        'checksum': "%s%s/md5sum.txt" % (ody_config.OUTPUT_DIR, wildcards.run),
        'fastqc': "%s%s/%s/fastqc.processed" % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR),
        'lims': "%s%s/%s/update_lims_db.processed" % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR),
        'sample_sheet': "%s%s/SampleSheet.csv" % (ody_config.OUTPUT_DIR, wildcards.run),
        'run_info': "%s%s/RunInfo.xml" % (ody_config.OUTPUT_DIR, wildcards.run)
    }
    input['demux'] = '%s%s/%s/demultiplex.processed' % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR)
    return input


rule publish:
    """
    copy all output to a published location
    """
    input:
        unpack(publish_input)
    output:
        touch(expand("{source}{{run}}/{status}/ody.complete", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    run:
        update_analysis({'step': 'publish', 'status': 'processing'})
        shell("rsync --info=STATS -rtl --perms --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r {ody_config.OUTPUT_DIR}{wildcards.run}/ {ody_config.PUBLISHED_DIR}{wildcards.run}/")
        send_success_email()

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    message = 'run %s failed\n see logs here: %s%s.log\n' % (config['run'], ody_config.LOG_DIR, config['run'])
    subject = 'Run Failed: %s' % config['run']
    sent = buildmessage(message, subject, {}, ody_config.EMAIL['from_email'], ody_config.EMAIL['admin_email'])
    update_analysis({'status': 'failed'})

def send_success_email():
    message = 'run %s completed successfully\n see logs here: %s%s.log\n' % (config['run'], ody_config.LOG_DIR, config['run'])
    cmd_file = '%s%s/script/demultiplex.sh' % (ody_config.OUTPUT_DIR, config['run'])
    cmd = util.get_file_contents(cmd_file)
    ss_file = '%s%s/SampleSheet.csv' % (ody_config.SOURCE_CLUSTER_PATH, config['run'])
    fastq_dir = '%s%s/fastq' % (ody_config.OUTPUT_CLUSTER_PATH, config['run'])
    summary_data = parse_stats.get_summary(fastq_dir, instrument, ss_file, config['run'])
    summary_data['cmd'] = cmd
    summary_data['version'] = 'bcl2fastq2 v2.2'
    subject = 'Demultiplex Summary for ' + config['run']
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL['from_email'], ody_config.EMAIL['to_email'], 'summary.html')
