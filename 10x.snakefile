'''
snakemake workflow for 10x single cell

Created on  2019-03-25

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2019 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

configfile: "/app/odybcl2fastq/snakemake_config.json"
localrules: all, update_lims_db, cp_source_to_output, checksum, publish, demultiplex_10x_cmd, count_10x_cmd, fastqc_cmd, fastq_email, insert_run_into_bauer_db
from odybcl2fastq.run import update_lims_db, setup_run_logger
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq import config as ody_config
from odybcl2fastq.bauer_db import BauerDB
import odybcl2fastq.util as util
import pandas as pd
import os

# parse sample sheet for sample names
sample_sheet_path = "%s%s/SampleSheet.csv" % (ody_config.SOURCE_DIR, config['run'])
with open(sample_sheet_path, 'r') as ln:
    idx = next(i for i, j in enumerate(ln) if j.startswith('[Data]'))
data = pd.read_csv(sample_sheet_path, skiprows=idx+1)

samples = list(data['Sample_ID'])
projects = list(data['Sample_Project'])
STATUS_DIR = 'status_test' if ody_config.TEST else 'status'

# set up bauer db for step updates
bauer = BauerDB(sample_sheet_path)

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

rule demultiplex_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("{source}{{run}}/{status}/analysis_id", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        run_dir=expand("{source}{{run}}/SampleSheet.csv", source=ody_config.SOURCE_DIR)
    output:
        expand("{output}{{run}}/script/demultiplex_10x.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        dual_index="--ignore-dual-index"
        if [ ! -z "{config[atac]}" ]; then
            dual_index=""
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq\n"
        cmd+="mkdir -p /scratch/{wildcards.run}_fastq_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.run}_fastq_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} mkfastq $dual_index --run={ody_config.SOURCE_CLUSTER_PATH}{wildcards.run} --samplesheet={ody_config.SOURCE_CLUSTER_PATH}{wildcards.run}/SampleSheet.csv --output-dir={ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq --localmem=\$((9*\$(ulimit -m)/10000000)) --loading-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --writing-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --processing-threads=\$SLURM_JOB_CPUS_PER_NODE --localcores=\$SLURM_JOB_CPUS_PER_NODE --barcode-mismatches=0 || exit_code=\$?\n"
        cmd+="cp -p */*.mri.tgz {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.run}_fastq_\$SLURM_JOB_ID\n"
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
        expand("{output}{{run}}/script/demultiplex_10x.sh", output=ody_config.OUTPUT_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

rule count_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        expand("{source}{{run}}/{status}/fastq_email.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR)
    output:
        expand("{output}{{run}}/script/{{project}}.{{sample}}_count.sh", output=ody_config.OUTPUT_DIR)
    shell:
        """
        fastq_path="{ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/fastq/{wildcards.project}/{wildcards.sample}"
        transcriptome="--transcriptome={ody_config.REF_PATH}{config[ref]}"
        if [ ! -z "{config[atac]}" ]; then
            transcriptome="--reference={ody_config.REF_PATH}{config[ref]}"
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/count/{wildcards.sample}\n"
        cmd+="mkdir -p /scratch/{wildcards.run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} count --id={wildcards.sample} $transcriptome --sample={wildcards.sample} --fastqs=$fastq_path --localmem=\$((9*\$(ulimit -m)/10000000)) --localcores=\$SLURM_JOB_CPUS_PER_NODE || exit_code=\$?\n\n"
        cmd+="/usr/bin/time -v cp -Rp {wildcards.sample}/*.mri.tgz {wildcards.sample}/outs {ody_config.OUTPUT_CLUSTER_PATH}{wildcards.run}/count/{wildcards.sample}/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.run}_{wildcards.sample}_\$SLURM_JOB_ID\n"
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
        script=expand("{output}{{run}}/script/{{project}}.{{sample}}_count.sh", output=ody_config.OUTPUT_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/{{project}}.{{sample}}_count.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    shell:
        """
        {input}
        """

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

rule fastq_email:
    """
    copy fastq files to final and send an email that fastq files are ready, this step is run only for jobs
    running count
    """
    input:
        checksum=expand("{output}{{run}}/md5sum.txt", output=ody_config.OUTPUT_DIR),
        fastqc=expand("{source}{{run}}/{status}/fastqc.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        lims=expand("{source}{{run}}/{status}/update_lims_db.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        demultiplex=expand("{source}{{run}}/{status}/demultiplex.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR),
        sample_sheet=expand("{output}{{run}}/SampleSheet.csv", output=ody_config.OUTPUT_DIR),
        run_info=expand("{output}{{run}}/RunInfo.xml", output=ody_config.OUTPUT_DIR)
    output:
        touch(expand("{source}{{run}}/{status}/fastq_email.processed", source=ody_config.SOURCE_DIR, status=STATUS_DIR))
    run:
        update_analysis({'step': 'count', 'status': 'processing'})
        shell("cp -r {ody_config.OUTPUT_DIR}{wildcards.run} {ody_config.PUBLISHED_DIR}")
        subject = 'Demultiplex Summary for: %s (count pending)' % config['run']
        send_success_email(subject)

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
    # if reference is provided then run count
    if config['ref']:
        # copy fastq to final and send an email that count is pending
        input['email'] = "%s%s/%s/fastq_email.processed" % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR)
        for i, sample in enumerate(samples):
            project = projects[i]
            key = 'count_%s.%s' % (project, sample)
            input[key] = "%s%s/%s/%s.%s_count.processed" % (ody_config.SOURCE_DIR, wildcards.run, STATUS_DIR, project, sample)
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
        subject = 'Demultiplex Summary for: %s' % config['run']
        send_success_email(subject)

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    message = 'run %s failed\n see logs here: %s%s.log\n' % (config['run'], ody_config.LOG_DIR, config['run'])
    subject = 'Run Failed: %s' % config['run']
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
    message = 'run %s completed successfully\n see logs here: %s%s.log\n' % (config['run'], ody_config.LOG_DIR, config['run'])
    cmd_file = '%s%s/script/demultiplex_10x.sh' % (ody_config.OUTPUT_DIR, config['run'])
    cmd = util.get_file_contents(cmd_file)
    ss_file = '%s%s/SampleSheet.csv' % (ody_config.SOURCE_DIR, config['run'])
    summary_data = get_summary_data(cmd, config['run'], ss_file)
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL['from_email'], ody_config.EMAIL['to_email'], '10x_summary.html')

def update_analysis(data):
    if not ody_config.TEST:
        analysis_file_path = '%s%s/%s/analysis_id' % (ody_config.SOURCE_DIR, config['run'], STATUS_DIR)
        if os.path.isfile(analysis_file_path):
            with open(analysis_file_path, 'r') as ln:
                analysis_id = ln.readline().strip()
            bauer.update_data('requests', analysis_id, data)
