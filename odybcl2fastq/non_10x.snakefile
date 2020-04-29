'''
snakemake workflow for non 10x

Created on  2020-03-26

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

include: "shared.snakefile"

localrules: all, update_lims_db, cp_source_to_output, checksum, publish, demultiplex_cmd, fastqc_cmd, insert_run_into_bauer_db
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq.parsers import parse_stats


MASK_SHORT_ADAPTER_READS = 22

# get some information about the run
sample_sheet = SampleSheet(sample_sheet_path)
instrument = sample_sheet.get_instrument()
run_info = '/source/' + config['run'] + '/RunInfo.xml'
run_type = sample_sheet.get_run_type()

# get basemasks for entire run and then use the one for the indexing strategy
# passed into config['mask_suffix'], for runs with only one strategy suffix will be
# an empty string
mask_lists, mask_samples = extract_basemasks(sample_sheet.sections['Data'], run_info, instrument, run_type)
jobs_tot = len(mask_lists)
if jobs_tot > 1:
    mask = config['mask_suffix'].replace('_', ',')
    sample_sheet_path = sample_sheet.write_new_sample_sheet(mask_samples[mask], config['mask_suffix'])
    sample_sheet = SampleSheet(sample_sheet_path)
    run_type = sample_sheet.get_run_type()
else:
    mask = next(iter(mask_lists))
mask_opt_list = mask_lists[mask]
mask_switch = '--use-bases-mask'
mask_opt = (mask_switch + ' ' + (' ' + mask_switch + ' ').join(mask_opt_list))

rule all:
    """
    final output of workflow
    """
    input:
        expand("/source/{run}/{status}/ody.complete", run=config['run'], status=status_dir)

def get_params_from_sample_sheet(sample_sheet):
    # users can add params to the end of the HEADER section of the sample sheet
    bcl_params = [
            '--no-lane-splitting',
            '--barcode-mismatches',
            '--mask-short-adapter-reads',
            '--minimum-trimmed-read-length',
            '--adapter-stringency',
            '--processing-threads',
            '--create-fastq-for-index-reads',
            '--mask-short-adapter-reads'
    ]
    params = {}
    for k, v in sample_sheet.sections['Header'].items():
        key = '--' + k.strip()
        if key in bcl_params:
            v = v.strip()
            params[key] = v if v else None
    return params

def shortest_read(r):
    return int(r[min(r.keys(), key=(lambda k:int(r[k])))])

def get_bcl_params(wildcards):
    # start with the defaults
    param_dict = {
        '--barcode-mismatches': 0,
        '--ignore-missing-bcls': None,
        '--ignore-missing-filter': None,
        '--ignore-missing-positions': None,
        '--minimum-trimmed-read-length': 1
    }
    if instrument in ['nextseq', 'miseq']:
        param_dict['--no-lane-splitting'] = None
    # check for short reads, do not mask
    if run_type == 'indrop' or shortest_read(sample_sheet.sections['Reads']) < MASK_SHORT_ADAPTER_READS:
        param_dict['--mask-short-adapter-reads'] = 0
    # grab any manually added params from sample sheet
    ss_params = get_params_from_sample_sheet(sample_sheet)
    param_dict.update(ss_params)
    param_list = [(k + ' ' + str(v)) if v is not None else k for k, v in param_dict.items()]
    param_str = ' '.join(param_list)
    return param_str

rule demultiplex_cmd:
    """
    build a bash file with the demux cmd
    """
    input:
        expand("/source/{run}/{status}/analysis_id", run=config['run'], status=status_dir),
        run_dir=expand("/source/{run}/SampleSheet{suffix}.csv", run=config['run'], suffix=config['suffix'])
    params:
        bcl_params=get_bcl_params
    output:
        expand("/analysis/{run}{suffix}/script/demultiplex.sh", run=config['run'], suffix=config['suffix'])
    shell:
        """
        cmd="#!/bin/bash\n"
        cmd+="ulimit -u \$(ulimit -Hu)\n"
        cmd+="exit_code=0\n"
        cmd+="mkdir -p {ody_config.OUTPUT_CLUSTER_PATH}{config[run]}{config[suffix]}/fastq\n"
        cmd+="/usr/bin/time -v bcl2fastq {params.bcl_params} --sample-sheet /source{config[run]}/SampleSheet{config[suffix]}.csv --runfolder-dir /source{config[run]} --output-dir {ody_config.OUTPUT_CLUSTER_PATH}{config[run]}{config[suffix]}/fastq --processing-threads 8 {mask_opt} || exit_code=\$?\n"
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
        expand("/analysis/{run}{suffix}/script/demultiplex.sh", run=config['run'], suffix=config['suffix'])
    output:
        touch(expand("/source/{run}/{status}/demultiplex.processed", run=config['run'], status=status_dir))
    run:
        update_analysis({'step': 'demultiplex', 'status': 'processing'})
        shell("{input}")

def publish_input(wildcards):
    """
    determine which files need to be ready to publish the run
    count is not run if there is not reference genome
    """
    input = {
        'demux': '/source/%s/%s/demultiplex.processed' % (config['run'], status_dir),
        'checksum': "/analysis/%s%s/md5sum.txt" % (config['run'], config['suffix']),
        'fastqc': "/source/%s/%s/fastqc.processed" % (config['run'], status_dir),
        'lims': "/source/%s/%s/update_lims_db.processed" % (config['run'], status_dir),
        'sample_sheet': "/analysis/%s%s/SampleSheet.csv" % (config['run'], config['suffix']),
        'run_info': "/analysis/%s%s/RunInfo.xml" % (config['run'], config['suffix'])
    }
    return input


rule publish:
    """
    copy all output to a published location
    """
    input:
        unpack(publish_input)
    output:
        touch(expand("/source/{{run}}/{status}/ody.complete", status=status_dir))
    run:
        update_analysis({'step': 'publish', 'status': 'processing'})
        shell("rsync --info=STATS -rtl --perms --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r /analysis/{config[run]}{config[suffix]}/ /published/{config[run]}{config[suffix]}/")
        send_success_email()

onsuccess:
    update_analysis({'status': 'complete'})

onerror:
    output_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s failed\n see logs here: %s%s.log\n' % (output_dir, ody_config.LOG_DIR, output_dir)
    subject = 'Run Failed: %s' % (output_dir)
    sent = buildmessage(message, subject, {}, ody_config.EMAIL['from_email'], ody_config.EMAIL['admin_email'])
    update_analysis({'status': 'failed'})

def send_success_email():
    output_dir = '%s%s' % (config['run'], config['suffix'])
    message = 'run %s completed successfully\n see logs here: %s%s.log\n' % (output_dir, ody_config.LOG_DIR, output_dir)
    cmd_file = '/analysis/%s/script/demultiplex.sh' % (output_dir)
    cmd = util.get_file_contents(cmd_file)
    ss_file = '/source%s/SampleSheet.csv' % (config['run'])
    fastq_dir = '%s%s/fastq' % (ody_config.OUTPUT_CLUSTER_PATH, output_dir)
    summary_data = parse_stats.get_summary(fastq_dir, instrument, ss_file, output_dir)
    summary_data['cmd'] = cmd
    summary_data['version'] = 'bcl2fastq2 v2.2'
    subject = 'Demultiplex Summary for ' + output_dir
    sent = buildmessage(message, subject, summary_data, ody_config.EMAIL['from_email'], ody_config.EMAIL['to_email'], 'summary.html')
