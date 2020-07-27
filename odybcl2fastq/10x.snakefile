'''
snakemake workflow for 10x single cell

Created on  2020-04-02

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@author: Nathan Weeks <nweeks@fas.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

include: "shared.snakefile"
localrules: all, update_lims_db, cp_source_to_output, checksum, publish, demultiplex_10x_cmd, count_10x_cmd, fastqc_cmd, multiqc, fastq_email, insert_run_into_bauer_db

sample_sheet = SampleSheet(sample_sheet_path)
samples = sample_sheet.get_samples()
projects = sample_sheet.get_projects()

rule all:
    """
    final output of workflow
    """
    input:
        expand("/data/source/{run}/{status}/ody.complete", run=config['run'], status=status_dir)

rule demultiplex_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("/data/source/{{run}}/{status}/analysis_id", status=status_dir),
        sample_sheet=expand("/data/source/{{run}}/SampleSheet{{suffix}}.csv")
    output:
        expand("/data/analysis/{{run}}{{suffix}}/script/demultiplex_10x.sh")
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p /data/analysis/{wildcards.run}{wildcards.suffix}/fastq\n"
        cmd+="mkdir -p /scratch/{wildcards.run}{wildcards.suffix}_fastq_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.run}{wildcards.suffix}_fastq_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} mkfastq --run=/data/source/{wildcards.run} --samplesheet=/data/source/{wildcards.run}/SampleSheet{wildcards.suffix}.csv --output-dir=/data/analysis/{wildcards.run}{wildcards.suffix}/fastq --localmem=\$((9*\$(ulimit -m)/10000000)) --loading-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --writing-threads=\$((SLURM_JOB_CPUS_PER_NODE/4)) --processing-threads=\$SLURM_JOB_CPUS_PER_NODE --localcores=\$SLURM_JOB_CPUS_PER_NODE --barcode-mismatches=0 || exit_code=\$?\n"
        cmd+="cp -p */*.mri.tgz /data/analysis/{wildcards.run}{wildcards.suffix}/fastq/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.run}{wildcards.suffix}_fastq_\$SLURM_JOB_ID\n"
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
        expand("/data/analysis/{{run}}{suffix}/script/demultiplex_10x.sh", suffix=config['suffix'])
    output:
        touch(expand("/data/source/{{run}}/{status}/demultiplex.processed", status=status_dir))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

rule count_10x_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("/data/source/{{run}}/{status}/demultiplex.processed", status=status_dir),
        expand("/data/source/{{run}}/{status}/fastq_email.processed", status=status_dir)
    output:
        expand("/data/analysis/{{run}}{{suffix}}/script/{{project}}.{{sample}}_count.sh")
    shell:
        """
        fastq_path="/data/analysis/{wildcards.run}{wildcards.suffix}/fastq/{wildcards.project}/{wildcards.sample}"
        transcriptome="--transcriptome=/ref/{config[ref]}"
        if [ ! -z "{config[atac]}" ]; then
            transcriptome="--reference=/ref/{config[ref]}"
        fi
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p /data/analysis/{wildcards.run}{wildcards.suffix}/count/{wildcards.sample}\n"
        cmd+="mkdir -p /scratch/{wildcards.run}{wildcards.suffix}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="cd /scratch/{wildcards.run}{wildcards.suffix}_{wildcards.sample}_\$SLURM_JOB_ID\n"
        cmd+="/usr/bin/time -v cellranger{config[atac]} count --id={wildcards.sample} $transcriptome --sample={wildcards.sample} --fastqs=$fastq_path --localmem=\$((9*\$(ulimit -m)/10000000)) --localcores=\$SLURM_JOB_CPUS_PER_NODE || exit_code=\$?\n\n"
        cmd+="/usr/bin/time -v cp -Rp {wildcards.sample}/*.mri.tgz {wildcards.sample}/outs /data/analysis/{wildcards.run}{wildcards.suffix}/count/{wildcards.sample}/ || exit_code=\$((exit_code | \$?))\n"
        cmd+="rm -rf /scratch/{wildcards.run}{wildcards.suffix}_{wildcards.sample}_\$SLURM_JOB_ID\n"
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
        script=expand("/data/analysis/{{run}}{suffix}/script/{{project}}.{{sample}}_count.sh", suffix=config['suffix'])
    output:
        touch(expand("/data/source/{{run}}/{status}/{{project}}.{{sample}}_count.processed", status=status_dir))
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
        checksum=expand("/data/analysis/{{run}}{suffix}/md5sum.txt", suffix=config['suffix']),
        fastqc=expand("/data/source/{{run}}/{status}/fastqc.processed", status=status_dir),
        multiqc=expand("/data/source/{{run}}/{status}/multiqc.processed", status=status_dir),
        lims=expand("/data/source/{{run}}/{status}/update_lims_db.processed", status=status_dir),
        demultiplex=expand("/data/source/{{run}}/{status}/demultiplex.processed", status=status_dir),
        sample_sheet=expand("/data/analysis/{{run}}{suffix}/SampleSheet.csv", suffix=config['suffix']),
        run_info=expand("/data/analysis/{{run}}{suffix}/RunInfo.xml", suffix=config['suffix'])
    output:
        touch(expand("/data/source/{{run}}/{status}/fastq_email.processed", status=status_dir))
    run:
        update_analysis({'step': 'count', 'status': 'processing'})
        # remove the published directory (if already exists) to avoid retaining any old files,
        shutil.rmtree(path="/data/published/{wildcards.run}{config[suffix]}", ignore_errors=True)
        # recursively hard-link analysis directory to published for speed & disk-usage reduction
        shutil.copytree(src="/data/analysis/{wildcards.run}{config[suffix]}",
                        dst="/data/published/{wildcards.run}{config[suffix]}",
                        symlinks=True, copy_function=os.link)
        shell("cp -r /data/analysis/{wildcards.run}{config[suffix]} /data/published/")
        subject = 'Demultiplex Summary for: %s%s (count pending)' % (wildcards.run, config['suffix'])
        send_success_email(subject)

def publish_input(wildcards):
    """
    determine which files need to be ready to publish the run
    count is not run if there is not reference genome
    """
    input = {
        'checksum': "/data/analysis/%s%s/md5sum.txt" % (wildcards.run, config['suffix']),
        'fastqc': "/data/source/%s/%s/fastqc.processed" % (wildcards.run, status_dir),
        'multiqc': "/data/source/%s/%s/multiqc.processed" % (wildcards.run, status_dir),
        'lims': "/data/source/%s/%s/update_lims_db.processed" % (wildcards.run, status_dir),
        'sample_sheet': "/data/analysis/%s%s/SampleSheet.csv" % (wildcards.run, config['suffix']),
        'run_info': "/data/analysis/%s%s/RunInfo.xml" % (wildcards.run, config['suffix'])
    }
    # if reference is provided then run count
    if config['ref']:
        # copy fastq to final and send an email that count is pending
        input['email'] = "/data/source/%s/%s/fastq_email.processed" % (wildcards.run, status_dir)
        for i, sample in enumerate(samples):
            project = projects[i]
            key = 'count_%s.%s' % (project, sample)
            input[key] = "/data/source/%s/%s/%s.%s_count.processed" % (wildcards.run, status_dir, project, sample)
    return input


rule publish:
    """
    copy all output to a published location
    """
    input:
        unpack(publish_input)
    output:
        touch(expand("/data/source/{{run}}/{status}/ody.complete", status=status_dir))
    run:
        update_analysis({'step': 'publish', 'status': 'processing'})
        # FIXME: need to use hard links...
        shell("rsync --info=STATS -rtl --perms --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r /data/analysis/{wildcards.run}{config[suffix]}/ /data/published/{wildcards.run}{config[suffix]}/")
        subject = 'Demultiplex Summary for: %s%s' % (wildcards.run, config['suffix'])
        send_success_email(subject)

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    run_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s failed\n see logs here: /data/log/%s.log\n' % (run_dir, run_dir)
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
    run_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s completed successfully\n see logs here: /data/log/%s.log\n' % (run_dir, run_dir)
    cmd_file = '/data/analysis/%s/script/demultiplex_10x.sh' % (run_dir)
    cmd = util.get_file_contents(cmd_file)
    summary_data = get_summary_data(cmd, run_dir, sample_sheet_path)
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL_FROM, ody_config.EMAIL_TO, '10x_summary.html')
