'''
snakemake workflow for 10x single cell

Created on  2020-04-02

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@author: Nathan Weeks <nweeks@fas.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

include: "shared.snakefile"
localrules: all, update_lims_db, cp_source_to_output, checksum, publish, demultiplex_10x_cmd, count_10x_cmd, fastqc_cmd, fastq_email, insert_run_into_bauer_db
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
import odybcl2fastq.util as util
import pandas as pd
import os

# parse sample sheet for sample names
with open(sample_sheet_path, 'r') as ln:
    idx = next(i for i, j in enumerate(ln) if j.startswith('[Data]'))
data = pd.read_csv(sample_sheet_path, skiprows=idx+1)

samples = list(data['Sample_ID'])
projects = list(data['Sample_Project'])

rule all:
    """
    final output of workflow
    """
    input:
        #expand("{source}{run}/{status}/ody.complete", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)
        #expand("{output}{run}{run_suffix}/script/demultiplex_10x.sh", output=ody_config.OUTPUT_DIR, run=config['run'], run_suffix=config['run_suffix'])
        expand("{source}{run}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)

rule demultiplex_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("{source}{{run}}/{status}/analysis_id", source=ody_config.SOURCE_DIR, status=status_dir),
        sample_sheet=expand("{source}{run}/SampleSheet{{run_suffix}}.csv", source=ody_config.SOURCE_DIR, run=config['run'])
    output:
        expand("{output}{{run}}{{run_suffix}}/script/demultiplex_10x.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        dual_index="--ignore-dual-index"
        if [ ! -z "{config[atac]}" ]; then
            dual_index=""
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}{wildcards.run_suffix}/fastq\n"
        cmd+="mkdir -p /scratch/{wildcards.run}{wildcards.run_suffix}_fastq_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.run}{wildcards.run_suffix}_fastq_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} mkfastq $dual_index --run={ody_config.SOURCE_CLUSTER_PATH}{wildcards.run} --samplesheet={ody_config.SOURCE_CLUSTER_PATH}{wildcards.run}/SampleSheet{wildcards.run_suffix}.csv --output-dir={ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}{wildcards.run_suffix}/fastq --localmem=\$((9*\$(ulimit -m)/10000000)) --loading-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --writing-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --processing-threads=\$SLURM_JOB_CPUS_PER_NODE --localcores=\$SLURM_JOB_CPUS_PER_NODE --barcode-mismatches=0 || exit_code=\$?\n"
        cmd+="cp -p */*.mri.tgz {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}{wildcards.run_suffix}/fastq/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.run}{wildcards.run_suffix}_fastq_\$SLURM_JOB_ID\n"
        cmd+="exit \$exit_code"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule demultiplex_10x:
    """
    run bash file for demux
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        expand("{output}{{run}}{run_suffix}/script/demultiplex_10x.sh", output=ody_config.OUTPUT_DIR, run_suffix=config['run_suffix'])
    output:
        touch(expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=status_dir))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

rule count_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("{source}{run}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir),
        expand("{source}{run}/{status}/fastq_email.processed", source=ody_config.SOURCE_DIR, run=config['run'], status=status_dir)
    output:
        expand("{output}{{output_run}}/script/{{project}}.{{sample}}_count.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        fastq_path="{ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/fastq/{wildcards.project}/{wildcards.sample}"
        transcriptome="--transcriptome={ody_config.REF_PATH}{config[ref]}"
        if [ ! -z "{config[atac]}" ]; then
            transcriptome="--reference={ody_config.REF_PATH}{config[ref]}"
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/count/{wildcards.sample}\n"
        cmd+="mkdir -p /scratch/{wildcards.output_run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.output_run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} count --id={wildcards.sample} $transcriptome --sample={wildcards.sample} --fastqs=$fastq_path --localmem=\$((9*\$(ulimit -m)/10000000)) --localcores=\$SLURM_JOB_CPUS_PER_NODE || exit_code=\$?\n\n"
        cmd+="/usr/bin/time -v cp -Rp {wildcards.sample}/*.mri.tgz {wildcards.sample}/outs {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.output_run}/count/{wildcards.sample}/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.output_run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="exit \$exit_code"
        echo "$cmd" >> {output}
        chmod 775 {output}
        """

rule count_10x:
    """
    run bash file for count
    the slurm_submit.py script will add slurm params to the top of this file
    """
    input:
        script=expand("{output}{output_run}/script/{{project}}.{{sample}}_count.sh", output=ody_config.OUTPUT_DIR, output_run=output_run)
    output:
        touch(expand("{source}{{run}}/{status}/{{project}}.{{sample}}_count.processed", source=ody_config.SOURCE_DIR, status=status_dir))
    shell:
        """
        {input}
        """

rule fastq_email:
    """
    copy fastq files to final and send an email that fastq files are ready, this step is run only for jobs
    running count
    """
    input:
        checksum=expand("{output}{output_run}/md5sum.txt", output=ody_config.OUTPUT_DIR, output_run=output_run),
        fastqc=expand("{source}{{run}}/{status}/fastqc.processed", source=ody_config.SOURCE_DIR, status=status_dir),
        lims=expand("{source}{{run}}/{status}/update_lims_db.processed", source=ody_config.SOURCE_DIR, status=status_dir),
        demultiplex=expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=status_dir),
        sample_sheet=expand("{output}{output_run}/{sample_sheet_name}", output=ody_config.OUTPUT_DIR, output_run=output_run, sample_sheet_name=sample_sheet_name),
        run_info=expand("{output}{output_run}/RunInfo.xml", output=ody_config.OUTPUT_DIR, output_run=output_run)
    output:
        touch(expand("{source}{{run}}/{status}/fastq_email.processed", source=ody_config.SOURCE_DIR, status=status_dir))
    run:
        update_analysis({'step': 'count', 'status': 'processing'})
        shell("cp -r {ody_config.OUTPUT_DIR}{wildcards.output_run} {ody_config.PUBLISHED_DIR}")
        subject = 'Demultiplex Summary for: %s (count pending)' % {wildcards.output_run}
        send_success_email(subject)

def publish_input(wildcards):
    """
    determine which files need to be ready to publish the run
    count is not run if there is not reference genome
    """
    input = {
        'checksum': "%s%s/md5sum.txt" % (ody_config.OUTPUT_DIR, output_run),
        'fastqc': "%s%s/%s/fastqc.processed" % (ody_config.SOURCE_DIR, wildcards.run, status_dir),
        'lims': "%s%s/%s/update_lims_db.processed" % (ody_config.SOURCE_DIR, wildcards.run, status_dir),
        'sample_sheet': "%s%s/SampleSheet.csv" % (ody_config.OUTPUT_DIR, output_run),
        'run_info': "%s%s/RunInfo.xml" % (ody_config.OUTPUT_DIR, wildcards.run)
    }
    # if reference is provided then run count
    if config['ref']:
        # copy fastq to final and send an email that count is pending
        input['email'] = "%s%s/%s/fastq_email.processed" % (ody_config.SOURCE_DIR, wildcards.run, status_dir)
        for i, sample in enumerate(samples):
            project = projects[i]
            key = 'count_%s.%s' % (project, sample)
            input[key] = "%s%s/%s/%s.%s_count.processed" % (ody_config.SOURCE_DIR, wildcards.run, status_dir, project, sample)
    return input


rule publish:
    """
    copy all output to a published location
    """
    input:
        unpack(publish_input)
    output:
        touch(expand("{source}{{run}}/{status}/ody.complete", source=ody_config.SOURCE_DIR, status=status_dir))
    run:
        update_analysis({'step': 'publish', 'status': 'processing'})
        shell("rsync --info=STATS -rtl --perms --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r {ody_config.OUTPUT_DIR}{output_run}/ {ody_config.PUBLISHED_DIR}{output_run}/")
        subject = 'Demultiplex Summary for: %s' % config['run']
        send_success_email(subject)

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    message = 'run %s failed\n see logs here: %s%s.log\n' % (output_run, ody_config.LOG_DIR, output_run)
    subject = 'Run Failed: %s' % output_run
    sent = buildmessage(message, subject, {}, ody_config.EMAIL['from_email'], ody_config.EMAIL['admin_email'])
    update_analysis({'status': 'failed'})

def get_summary_data(cmd, run, ss_file):
    sample_sheet = util.get_file_contents(ss_file)
    assumptions = [
        'chromium single-cell RNA-seq',
        'ignoring the second index on a dual-indexed flowcell'
    ]
    versions = [
        'cellranger 3.1.0',
        'fastqc 0.11.8'
    ]
    if config['ref']:
        assumptions.append('reference genome %s, see annotation under versions below' % os.path.basename(config['ref']))
        assumptions.append('samples are full cell, not nuclei')
        versions.append('annotation gtf: %s' % config['gtf'])
    summary_data = {
        'fastq_url': ody_config.FASTQ_URL,
        'fastq_dir': ody_config.PUBLISHED_CLUSTER_PATH,
        'versions': versions,
        'run': run,
        'run_folder': run,
        'cmd': cmd,
        'sample_sheet_file': ss_file,
        'sample_sheet': sample_sheet,
        'assumptions': assumptions
    }
    return summary_data

def send_success_email(subject):
    message = 'run %s completed successfully\n see logs here: %s%s.log\n' % (output_run, ody_config.LOG_DIR, output_run)
    cmd_file = '%s%s/script/demultiplex_10x.sh' % (ody_config.OUTPUT_DIR, output_run)
    cmd = util.get_file_contents(cmd_file)
    summary_data = get_summary_data(cmd, output_run, sample_sheet_path)
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL['from_email'], ody_config.EMAIL['to_email'], '10x_summary.html')