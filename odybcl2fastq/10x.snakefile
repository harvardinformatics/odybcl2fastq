'''
snakemake workflow for 10x single cell

Created on  2020-04-02

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@author: Nathan Weeks <nweeks@fas.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

include: "shared.snakefile"
localrules: all, update_lims_db, cp_source_to_output, checksum_cmd, publish, demultiplex_10x_cmd, count_10x_cmd, fastqc_cmd, multiqc, fastq_email, insert_run_into_bauer_db

sample_sheet = SampleSheet(sample_sheet_path)
samples = sample_sheet.get_samples()
projects = sample_sheet.get_projects()

rule all:
    """
    final output of workflow
    """
    input:
        expand("/sequencing/source/{run}/{status}/ody.complete", run=config['run'], status=status_dir)

rule demultiplex_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        f"/sequencing/source/{config['run']}/{status_dir}/analysis_id",
        sample_sheet_path
    output:
        expand("/sequencing/analysis/{{run}}{{suffix}}/script/demultiplex_10x.sh")
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -n \$(ulimit -Hn)\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p /sequencing/analysis/{config[run]}{config[suffix]}/fastq\n"
        cmd+="mkdir -p /scratch/{config[run]}{config[suffix]}_fastq_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{config[run]}{config[suffix]}_fastq_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} mkfastq --run=/sequencing/source/{config[run]} --samplesheet={sample_sheet_path} --output-dir=/sequencing/analysis/{config[run]}{config[suffix]}/fastq --localmem=\$((9*\$(ulimit -m)/10000000)) --loading-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --writing-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --processing-threads=\$SLURM_JOB_CPUS_PER_NODE --localcores=\$SLURM_JOB_CPUS_PER_NODE --barcode-mismatches=0 || exit_code=\$?\n"
        cmd+="cp -p */*.mri.tgz /sequencing/analysis/{config[run]}{config[suffix]}/fastq/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{config[run]}{config[suffix]}_fastq_\$SLURM_JOB_ID\n"
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
        expand("/sequencing/analysis/{{run}}{suffix}/script/demultiplex_10x.sh", suffix=config['suffix'])
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/demultiplex.processed", status=status_dir))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

rule count_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("/sequencing/source/{{run}}/{status}/demultiplex.processed", status=status_dir),
        expand("/sequencing/source/{{run}}/{status}/fastq_email.processed", status=status_dir)
    output:
        expand("/sequencing/analysis/{{run}}{{suffix}}/script/{{project}}.{{sample}}_count.sh")
    shell:
        """
        fastq_path="/sequencing/analysis/{config[run]}{config[suffix]}/fastq"
        transcriptome="--transcriptome=/ref/{config[ref]}"
        if [ ! -z "{config[atac]}" ]; then
            transcriptome="--reference=/ref/{config[ref]}"
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p /sequencing/analysis/{config[run]}{config[suffix]}/count/{wildcards.sample}\n"
        cmd+="mkdir -p /scratch/{config[run]}{config[suffix]}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{config[run]}{config[suffix]}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} count --id={wildcards.sample} $transcriptome --sample={wildcards.sample} --fastqs=$fastq_path --localmem=\$((9*\$(ulimit -m)/10000000)) --localcores=\$SLURM_JOB_CPUS_PER_NODE || exit_code=\$?\n\n"
        cmd+="/usr/bin/time -v cp -Rp {wildcards.sample}/*.mri.tgz {wildcards.sample}/outs /sequencing/analysis/{config[run]}{config[suffix]}/count/{wildcards.sample}/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{config[run]}{config[suffix]}_{wildcards.sample}_\$SLURM_JOB_ID\n"
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
        script=expand("/sequencing/analysis/{{run}}{suffix}/script/{{project}}.{{sample}}_count.sh", suffix=config['suffix'])
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/{{project}}.{{sample}}_count.processed", status=status_dir))
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
        checksum=expand("/sequencing/analysis/{{run}}{suffix}/md5sum.txt", suffix=config['suffix']),
        fastqc=expand("/sequencing/source/{{run}}/{status}/fastqc.processed", status=status_dir),
        multiqc=expand("/sequencing/source/{{run}}/{status}/multiqc.processed", status=status_dir),
        lims=expand("/sequencing/source/{{run}}/{status}/update_lims_db.processed", status=status_dir),
        demultiplex=expand("/sequencing/source/{{run}}/{status}/demultiplex.processed", status=status_dir),
        sample_sheet=expand("/sequencing/analysis/{{run}}{suffix}/SampleSheet.csv", suffix=config['suffix']),
        run_info=expand("/sequencing/analysis/{{run}}{suffix}/RunInfo.xml", suffix=config['suffix'])
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/fastq_email.processed", status=status_dir))
    run:
        update_analysis({'step': 'count', 'status': 'processing'})
        # remove the published directory (if already exists) to avoid retaining any old files,
        shutil.rmtree(path=f"/sequencing/published/{config['run']}{config['suffix']}", ignore_errors=True)
        # recursively hard-link analysis directory to published for speed & disk-usage reduction
        shutil.copytree(src=f"/sequencing/analysis/{config['run']}{config['suffix']}",
                        dst=f"/sequencing/published/{config['run']}{config['suffix']}",
                        symlinks=True, copy_function=util.link_readonly)
        subject = 'Demultiplex Summary for: %s%s (count pending)' % (config['run'], config['suffix'])
        send_success_email(subject)

def publish_input(wildcards):
    """
    determine which files need to be ready to publish the run
    count is not run if there is not reference genome
    """
    input = {
        'checksum': "/sequencing/analysis/%s%s/md5sum.txt" % (config['run'], config['suffix']),
        'fastqc': "/sequencing/source/%s/%s/fastqc.processed" % (config['run'], status_dir),
        'multiqc': "/sequencing/source/%s/%s/multiqc.processed" % (config['run'], status_dir),
        'lims': "/sequencing/source/%s/%s/update_lims_db.processed" % (config['run'], status_dir),
        'sample_sheet': "/sequencing/analysis/%s%s/SampleSheet.csv" % (config['run'], config['suffix']),
        'run_info': "/sequencing/analysis/%s%s/RunInfo.xml" % (config['run'], config['suffix'])
    }
    # if reference is provided then run count
    if config['ref']:
        # copy fastq to final and send an email that count is pending
        input['email'] = "/sequencing/source/%s/%s/fastq_email.processed" % (config['run'], status_dir)
        for i, sample in enumerate(samples):
            project = projects[i]
            key = 'count_%s.%s' % (project, sample)
            input[key] = "/sequencing/source/%s/%s/%s.%s_count.processed" % (config['run'], status_dir, project, sample)
    return input


rule publish:
    """
    copy all output to a published location
    """
    input:
        unpack(publish_input)
    output:
        touch(expand("/sequencing/source/{{run}}/{status}/ody.complete", status=status_dir))
    run:
        update_analysis({'step': 'publish', 'status': 'processing'})
        shutil.copytree(src=f"/sequencing/analysis/{config['run']}{config['suffix']}",
                        dst=f"/sequencing/published/{config['run']}{config['suffix']}",
                        symlinks=True, copy_function=util.link_readonly, dirs_exist_ok=True)
        subject = 'Demultiplex Summary for: %s%s' % (config['run'], config['suffix'])
        send_success_email(subject)

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    run_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s failed\n see logs here: /sequencing/log/%s.log\n' % (run_dir, run_dir)
    subject = 'Run Failed: %s' % run_dir
    sent = buildmessage(message, subject, {}, ody_config.EMAIL_FROM, ody_config.EMAIL_TO)
    update_analysis({'status': 'failed'})

def get_summary_data(cmd, run, ss_file):
    sample_sheet = util.get_file_contents(ss_file)
    assumptions = [
        'chromium single-cell RNA-seq'
    ]
    versions = [
        'cellranger 4.0.0',
        'fastqc 0.11.8'
        'multiqc 1.9'
    ]
    if config['ref']:
        assumptions.append('reference genome %s, see annotation under versions below' % os.path.basename(config['ref']))
        versions.append('annotation gtf: %s' % config['gtf'])
    summary_data = {
        'fastq_url': ody_config.FASTQ_URL,
        'globus_url': ody_config.GLOBUS_URL,
        'fastq_dir': ody_config.PUBLISHED_CLUSTER_PATH,
        'versions': versions,
        'run': run,
        'cmd': cmd,
        'sample_sheet_file': ss_file,
        'sample_sheet': sample_sheet,
        'assumptions': assumptions
    }
    return summary_data

def send_success_email(subject):
    run_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s completed successfully\n see logs here: /sequencing/log/%s.log\n' % (run_dir, run_dir)
    cmd_file = '/sequencing/analysis/%s/script/demultiplex_10x.sh' % (run_dir)
    cmd = util.get_file_contents(cmd_file)
    summary_data = get_summary_data(cmd, run_dir, sample_sheet_path)
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL_FROM, ody_config.EMAIL_TO, '10x_summary.html')
