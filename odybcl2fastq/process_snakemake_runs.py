#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and process them

Created on  2020-04-02

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time, sys
import logging
import subprocess
import json
import traceback
from time import sleep
from datetime import datetime
from multiprocessing import Pool
from subprocess import Popen, PIPE, STDOUT
from odybcl2fastq import config, initLogger, setupMainLogger
from odybcl2fastq import constants as const
import odybcl2fastq.util as util
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from snakemake import snakemake
from odybcl2fastq.status_db import StatusDB

STATUS_DIR = 'status_test' if config.TEST else 'status'
PROCESSED_FILE_NAME = 'ody.processed'
PROCESSED_FILE = '%s/%s' % (STATUS_DIR, PROCESSED_FILE_NAME)
COMPLETE_FILE_NAME = 'ody.complete'
COMPLETE_FILE = '%s/%s' % (STATUS_DIR, COMPLETE_FILE_NAME)
SKIP_FILE = 'odybcl2fastq.skip'
INCOMPLETE_NOTIFIED_FILE = '%s/ody.incomplete_notified' % STATUS_DIR
DAYS_TO_SEARCH = 3
INCOMPLETE_AFTER_DAYS = 4
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('Aug 10 2019', '%b %d %Y')
REQUIRED_FILES = ['InterOp/QMetricsOut.bin', 'InterOp/TileMetricsOut.bin', 'RunInfo.xml', 'RTAComplete.txt']
TYPES_10X = ['10x single cell', '10x single cell rna', '10x single nuclei rna', '10x single cell vdj', '10x single cell atac']
PROC_NUM = int(os.getenv('ODYBCL2FASTQ_PROC_NUM', 7))
FREQUENCY = 60

logger = setupMainLogger()

def failure_email(run, log, cmd, ret_code, std_out, std_err = ''):
    subject = "Run Failed: %s" % run
    message = (
        "%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
        " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log)
    )
    send_email(message, subject)

def send_email(message, subject, to = 'to_email_error'):
    fromaddr = config.EMAIL_FROM
    toemaillist = config.EMAIL_TO
    buildmessage(message, subject, None, fromaddr, toemaillist)

def need_to_process(dir):
    now = datetime.now()
    m_time = datetime.fromtimestamp(os.stat(dir).st_mtime)
    # filter out if modified before cutover to odybcl2fastq
    if m_time < SEARCH_AFTER_DATE:
        return False
    # filter out if modified outside or search window
    if ((now - m_time).days) > DAYS_TO_SEARCH:
        return False
    # filter out if tagged as processed
    if os.path.isfile(dir + PROCESSED_FILE):
        return False
    # filter out if tagged as skip
    if os.path.isfile(dir + SKIP_FILE):
        return False
    # filter out if any required files are missing
    for req in REQUIRED_FILES:
        if not os.path.exists(dir + req):
            return False
    return True

def run_is_incomplete(dir):
    now = datetime.now()
    m_time = datetime.fromtimestamp(os.stat(dir).st_mtime)
    # filter out if modified before cutover to odybcl2fastq
    if m_time < SEARCH_AFTER_DATE:
        return False
    # filter out if modified after reasonable delay to allow for completion
    if ((now - m_time).days) <= INCOMPLETE_AFTER_DAYS:
        return False
    # filter out if tagged as complete
    if os.path.isfile(dir + COMPLETE_FILE):
        return False
    # filter out if never tagged for processing
    if not os.path.isfile(dir + PROCESSED_FILE):
        return False
    # filter out already notified
    if os.path.isfile(dir + INCOMPLETE_NOTIFIED_FILE):
        return False
    return True

def find_runs(filter):
    # get all subdirectories
    dirs = sorted(glob.glob('/source/*/'))
    runs = []
    for dir in dirs:
        if filter(dir):
            runs.append(dir)
    return runs

def check_sample_sheet(sample_sheet, run):
    # if sample sheet is not already there then copy the one from run_folder named
    # for flowcell
    if not os.path.exists(sample_sheet):
        flowcell = run.split('_')[-1][1:]
        path = '/source/sample_sheet/' + flowcell + '.csv'
        if os.path.exists(path):
            util.copy(path, sample_sheet)

def get_sample_sheet_path(run_dir):
    # set default
    sample_sheet_path = run_dir + 'SampleSheet.csv'
    # see if a txt file indicates a specific, existing sample sheet
    sample_sheet_txt = glob.glob(run_dir + 'SampleSheet*txt')
    # if there are more than one then just use default
    if len(sample_sheet_txt) == 1:
        sample_sheet_path_tmp = sample_sheet_txt[0].replace('.txt', '.csv')
        # override with this path if it exists
        if os.path.exists(sample_sheet_path_tmp):
            sample_sheet_path = sample_sheet_path_tmp
    return sample_sheet_path

def check_complete(run_dir):
    if not os.path.exists('%s%s' % (run_dir, COMPLETE_FILE)):
        sub_dirs = next(os.walk('%s%s' % (run_dir, STATUS_DIR)))[1]
        complete = True
        for sub_dir in sub_dirs:
            if not os.path.exists('%s%s/%s/%s' % (run_dir, STATUS_DIR, sub_dir, COMPLETE_FILE_NAME)):
                complete = False
        if complete:
            util.touch(run_dir, COMPLETE_FILE)

def get_custom_suffix(sample_sheet_path):
    suffix = ''
    if 'SampleSheet_' in sample_sheet_path:
        suffix = os.path.basename(sample_sheet_path).replace('SampleSheet_', '').replace('.csv', '')
    return suffix

def get_reference(run_dir, run_type, sample_sheet):
    run = os.path.basename(os.path.normpath(run_dir))
    subs = sample_sheet.get_submissions()
    stdb = StatusDB()
    # TODO: assuming that all samples in the run are same type if not we
    # should reconsider how to organize all the count bits
    ref = ''
    if subs:
        sams = stdb.minilims_select('Sample', None, 'Submission', subs[0])
        if sams:
            ref = stdb.minilims_select('Sample', sams[0][0], 'Reference_Genome')
            if ref:
                ref = ref[0][3]
    ref_file = ''
    gtf = ''
    if ref == 'hg19' or ref == 'human_hg19': # human
        ref_file = 'refdata-cellranger-hg19-3.0.0'
    elif ref == 'GRCh' or ref == 'human_GRC38': # human
        if run_type == '10x single cell atac':
            ref_file = 'atac-seq/refdata-cellranger-atac-GRCh38-1.0.1'
        else:
            ref_file = 'refdata-cellranger-GRCh38-3.0.0'
    elif 'Zebrafish' in ref or ref =='Zebrafish_GRCz11':
        if run_type == '10x single cell atac':
            ref_file = 'atac-seq/refdata-cellranger-atac-zebrafish/Danio_rerio.GRCz11'
        else:
            ref_file = 'zebrafish_ensembl/Danio_rerio.GRCz11'
    elif 'mouse' in ref or ref == 'mouse_mm10':
        if run_type == '10x single cell atac':
            ref_file = 'atac-seq/refdata-cellranger-atac-mm10-1.0.1'
        else:
            ref_file = 'refdata-cellranger-mm10-3.0.0'
    elif ref not in ['', 'None', 'Other']: # this will cause count to error and then we can add a genome
        # email admins to notify we need a reference genome
        message = "run %s doesn't have a reference genome prepared for: %s\n" % (run, ref)
        subject = 'Run needs reference genome: %s' % run
        sent = buildmessage(message, subject, {}, config.EMAIL_FROM, config.EMAIL_ADMIN)
    if ref_file:
        # get gtf file
        with open(('/ref/%s/reference.json' % (ref_file)), 'r') as f:
            data = json.load(f)
            if 'input_gtf_files' in data:
                gtf = ', '.join(data['input_gtf_files']).replace('.filtered.gtf', '').replace('.gtf.filtered', '')
    return (ref_file, gtf)

def get_run_suffix(custom_suffix, mask_suffix):
    suffix = ''
    if mask_suffix:
        suffix += '_%s' % mask_suffix
    elif custom_suffix:
        suffix += '_%s' % custom_suffix
    return suffix

def get_output_log(run):
    logdir = os.environ.get('ODYBCL2FASTQ_RUN_LOG_DIR', '/log/')
    return os.path.join(logdir, run + '.log')

def get_10x_snakemake_config(run_dir, run_type, sample_sheet, run, suffix):
    atac = ''
    if 'atac' in run_type:
        atac = '-atac'
    ref_file = ''
    gtf = ''
    if not 'nuclei' in run_type and not run_type == '10x single cell vdj':
        ref_file, gtf = get_reference(run_dir, run_type, sample_sheet)
    return {'run': run, 'ref': ref_file, 'gtf': gtf, 'atac': atac, 'suffix': suffix}

def get_ody_snakemake_opts(run_dir, ss_path, run_type, suffix, mask_suffix):
    run = os.path.basename(os.path.normpath(run_dir))
    sample_sheet = SampleSheet(ss_path)
    sample_sheet.validate()
    if run_type in TYPES_10X:
        snakemake_config = get_10x_snakemake_config(run_dir, run_type, sample_sheet, run, suffix)
        snakefile = '10x.snakefile'
    else:
        snakemake_config = {'run': run, 'suffix': suffix, 'mask_suffix': mask_suffix}
        snakefile = 'non_10x.snakefile'

    opts = {
        '--cores': 99,
        '--local-cores': 4,
        '--max-jobs-per-second': 2,
        '--config': ' '.join(['%s=%s' % (k, v) for k, v in snakemake_config.items()]),
        '--cluster-config': '/app/odybcl2fastq/snakemake_cluster.json',
        '--cluster': "'python /app/odybcl2fastq/slurm_submit.py'",
        '--cluster-status': "'python /app/odybcl2fastq/cluster_status.py'",
        '--printshellcmds': None,
        '--reason': None,
        '-s': '/app/odybcl2fastq/%s' % snakefile,
        '--jobscript': '/app/odybcl2fastq/jobscript.sh',
        '-d': '/snakemake/ody',
        '--latency-wait': 120
    }
    return [k + ((' %s' % v) if v else '') for k, v in opts.items()]

def run_snakemake(cmd, output_log):
    # run unix cmd, stream out and error, return last lines of out
    with open(output_log, 'a+') as writer:
        with open(output_log, 'r') as reader:
            proc = Popen(cmd, shell=True, stderr=STDOUT, stdout=writer)
            lines = []
            # stream output to log, std out
            while proc.poll() is None:
                line = reader.readline()
                # save last 40 lines for email
                if line:
                    lines.append(line)
                    if len(lines) > 40:
                        lines.pop(0)
            # write any remaining output
            try:
                line = reader.readline()
                lines.append(line)
            except Exception:
                pass
            code = proc.wait()
    return code, '\n'.join(lines)

def notify_incomplete_runs():
    run_dirs = find_runs(run_is_incomplete)
    run_dirs_str = "\n".join(run_dirs)
    if run_dirs:
        message = "The following runs failed to complete %s or more days ago:\n\n%s" % (INCOMPLETE_AFTER_DAYS, run_dirs_str)
        send_email(message, 'Odybcl2fastq incomplete runs')
        for run in run_dirs:
            util.touch(run, INCOMPLETE_NOTIFIED_FILE)


def tail(f, n):
    cmd = "tail -n %i %s | grep returned" % (n, f)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    std_out, std_err = proc.communicate()
    return std_out.splitlines(True)

def get_runs():
    run_dirs_tmp = find_runs(need_to_process)
    run_dirs = []
    for run_dir in run_dirs_tmp:
        run_info = {}
        try:
            # get the custom suffix so we only look at the sample sheet that is
            # flagged if one exists
            ss_path = get_sample_sheet_path(run_dir)
            custom_suffix = get_custom_suffix(ss_path)
            ss_path = get_sample_sheet_path('%s%s' % (run_dir, custom_suffix))
            run = os.path.basename(os.path.normpath(run_dir))
            # copy samplesheet from samplesheet folder if necessary
            check_sample_sheet(ss_path, run)
            sample_sheet = SampleSheet(ss_path)
            instrument = sample_sheet.get_instrument()
            sam_types = sample_sheet.get_sample_types()
            run_info_file = run_dir + '/RunInfo.xml'
            # create a list of runs with their type
            for t, v in sam_types.items():
                mask_suffix = ''
                if v not in TYPES_10X: #non 10x
                    # get some information about the run
                    # consider if multiple indexing strategies are needed
                    # meaning it will be run more than once
                    mask_lists, mask_samples = extract_basemasks(sample_sheet.sections['Data'], run_info_file, instrument, v, False)
                    jobs_tot = len(mask_lists)
                    if jobs_tot > 1:
                        for mask, mask_list in mask_lists.items():
                            mask_suffix = mask.replace(',', '_')
                            mask_status_path = '%s%s/%s/' % (run_dir, STATUS_DIR, mask_suffix)
                            # only start the run if the processed file for that mask
                            # is not present, this will enable restart of one mask
                            if not os.path.isfile(mask_status_path + PROCESSED_FILE_NAME):
                                run_dirs.append({'run':run_dir, 'type':v, 'mask_suffix': mask_suffix, 'custom_suffix': custom_suffix})
                                util.touch(mask_status_path, PROCESSED_FILE_NAME)
                    else:
                        run_dirs.append({'run':run_dir, 'type':v, 'mask_suffix': mask_suffix, 'custom_suffix': custom_suffix})
                else:
                    run_dirs.append({'run':run_dir, 'type':v, 'mask_suffix': mask_suffix, 'custom_suffix': custom_suffix})
                break
        except:
            traceback.print_exc()
    return run_dirs

def process_runs(pool):
    '''
    Fills up the pool with runs with apply_async
    If anything is ready, result is placed in success_runs or failed_runs
    Then looks for more runs
    '''
    logger.info("Processing runs")
    run_dirs = get_runs()
    logger.info("Found %s runs: %s\n" % (len(run_dirs), json.dumps(run_dirs)))

    results = {}
    failed_runs = []
    success_runs = []
    queued_runs = {}
    while len(run_dirs) > 0 or len(results) > 0:
        logger.info("Current runs found %s Runs pending results: %s Runs queued: %s\n" % (json.dumps(run_dirs), json.dumps(list(results.keys())),
                    json.dumps(list(queued_runs.keys()))))
        while len(run_dirs) > 0:
            run_info = run_dirs.pop()
            run_type = run_info['type']
            run_dir = run_info['run']
            mask_suffix = run_info['mask_suffix']
            custom_suffix = run_info['custom_suffix']
            run = '%s%s%s' % (os.path.basename(os.path.normpath(run_dir)), ('_' if mask_suffix else ''), mask_suffix)
            ss_path = get_sample_sheet_path(run_dir)
            suffix = get_run_suffix(custom_suffix, mask_suffix)
            opts = get_ody_snakemake_opts(run_dir, ss_path, run_type, suffix, mask_suffix)
            logger.info("Queueing odybcl2fastq cmd for %s:\n" % (run))
            run_log = get_output_log(run)
            cmd = 'snakemake ' + ' '.join(opts)
            msg = "Running cmd: %s\n" % cmd
            logger.info(msg)
            runlogger = initLogger('run_logger', run + '.log')
            runlogger.info(msg)
            # create status dir if it doesn't exist
            status_path = '%s%s' % (run_dir, STATUS_DIR)
            if not os.path.exists(status_path):
                os.mkdir(status_path)
            # touch processed file so errors in snakemake don't cause the run to
            # continually be queued
            util.touch(run_dir, PROCESSED_FILE)
            results[run] = pool.apply_async(run_snakemake, (cmd, run_log,))
            queued_runs[run] = 1
        for run in list(results.keys()):
            result = results[run]
            if result.ready():
                # 5 day timeout
                timeout = (5 * 24 * 60 * 60)
                ret_code, lines = result.get(timeout)
                if ret_code == 0:
                    success_runs.append(run)
                    status = 'success'
                    check_complete(run_dir)
                else:
                    failed_runs.append(run)
                    status = 'failure'
                    failure_email(run, run_log, cmd, ret_code, lines)
                    logging.info('Run failed: %s with code %s\n %s\n %s' % (run, str(ret_code), cmd, lines))
                del results[run]
                del queued_runs[run]

        logger.info("After checking results, runs found %s Runs pending results: %s Runs queued: %s\n"
                % (json.dumps(run_dirs), json.dumps(list(results.keys())),
                    json.dumps(list(queued_runs.keys()))))
        sleep(10)
        new_run_dirs = get_runs()
        for new_run_dir in new_run_dirs:
            new_run = os.path.basename(os.path.normpath(new_run_dir['run']))
            if new_run not in queued_runs:
                run_dirs.append(new_run_dir)
        if len(run_dirs) > 0:
            logger.info("Found %s more runs: %s\n" % (len(run_dirs),
                json.dumps(run_dirs)))
    logger.info(
        "Completed %s runs: %s success %s and %s failures %s\n\n\n" %
        ((len(success_runs) + len(failed_runs)), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs))
    )

def main():
    try:
        logger.info("Starting ody10x processing")
        logger.info("Running with ")
        proc_num = PROC_NUM
        pool = Pool(proc_num)
        # run continuously
        while True:
            # queue new runs for demultiplexing with bcl2fastq2
            process_runs(pool)
            # check for any runs that started but never completed demultiplexing
            notify_incomplete_runs()
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logger.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
    except Exception as e:
        logging.exception(e)
        #send_email(str(e), 'Odybcl2fastq exception')
        return 1


if __name__ == "__main__":
    sys.exit(main())
