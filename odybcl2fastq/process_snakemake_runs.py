#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina 10x runs and process them

Created on  2019-03-25

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2019 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time, sys
import logging
import subprocess
import json
from time import sleep
from datetime import datetime
from multiprocessing import Pool
from odybcl2fastq import config
from odybcl2fastq import constants as const
import odybcl2fastq.util as util
from odybcl2fastq import run as ody_run
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.parsers.samplesheet import SampleSheet
from snakemake import snakemake
from odybcl2fastq.status_db import StatusDB

LOG_HTML = config.FINAL_DIR + 'odybcl2fastq_log.html'
PROCESSED_FILE = 'status/ody.processed'
COMPLETE_FILE = 'status/ody.complete'
SKIP_FILE = 'odybcl2fastq.skip'
INCOMPLETE_NOTIFIED_FILE = 'odybcl2fastq.incomplete_notified'
DAYS_TO_SEARCH = 3
INCOMPLETE_AFTER_DAYS = 2
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('Mar 20 2019', '%b %d %Y')
REQUIRED_FILES = ['InterOp/QMetricsOut.bin', 'InterOp/TileMetricsOut.bin', 'RunInfo.xml', 'RTAComplete.txt']
PROC_NUM = int(os.getenv('ODYBCL2FASTQ_PROC_NUM', 2))

FREQUENCY = 3

logger = logging.getLogger('odybcl2fastq10x')

def get_logger():
    return logger

def failure_email(run, cmd, ret_code, std_out, std_err):
    log = ody_run.get_output_log(run)
    subject = "Run Failed: %s" % run
    message = (
        "%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
        " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log)
    )
    send_email(message, subject)


def send_email(message, subject):
    fromaddr = config.EMAIL['from_email']
    toemaillist = config.EMAIL['to_email']
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
    dirs = sorted(glob.glob(config.SOURCE_DIR + '*/'))
    runs = []
    for dir in dirs:
        if filter(dir):
            runs.append(dir)
    return runs


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

def get_reference(run_dir):
    ss_path = get_sample_sheet_path(run_dir)
    run = os.path.basename(os.path.normpath(run_dir))
    sample_sheet = SampleSheet(ss_path)
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
    if ref == 'hg19': # human
        ref_file = '%srefdata-cellranger-hg19-3.3.0' % config.ody['ref_dir']
    elif ref == 'GRCh': # human
        ref_file = '%srefdata-cellranger-GRCh38-3.0.0' % config.ody['ref_dir']
    elif 'Zebrafish' in ref:
        ref_file = '%szebrafish_ensembl' % (config.ody['ref_dir'])
    elif 'mouse' in ref:
        ref_file = '%srefdata-cellranger-mm10-3.0.0' % (config.ody['ref_dir'])
    else: # this will cause count to error and then we can add a genome
        # email admins to notify we need a reference genome
        message = "run %s doesn't have a reference genome prepared for: %s\n" % (run, ref)
        subject = 'Run needs reference genome: %s' % run
        sent = buildmessage(message, subject, {}, config.EMAIL['from_email'], config.EMAIL['admin_email'])
    # get gtf file
    with open(ref_file + '/reference.json', 'r') as f:
        data = json.load(f)
        if 'input_gtf_files' in data:
            gtf = ', '.join(data['input_gtf_files']).replace('.filtered.gtf', '').replace('.gtf.filtered', '')
    return (ref_file, gtf)

def get_output_log(run):
    logdir = os.environ.get('ODYBCL2FASTQ_RUN_LOG_DIR', config.LOG_DIR)
    return os.path.join(logdir, run + '.log')

def setup_run_logger(run):
    # take level from env or INFO
    runlogger = logging.getLogger('run_logger')
    level = logging.getLevelName(os.environ.get('ODYBCL2FASTQ_RUN_LOG_LEVEL', 'INFO'))
    runlogger.setLevel(level)
    handler = logging.FileHandler(get_output_log(run))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    handler.setLevel(level)
    handler.setFormatter(formatter)
    runlogger.addHandler(handler)
    return runlogger

def get_ody_snakemake_opts(run_dir):
    run = os.path.basename(os.path.normpath(run_dir))
    ss_path = get_sample_sheet_path(run_dir)
    sample_sheet = SampleSheet(ss_path)
    df = sample_sheet.get_samples()
    samples = df['Sample_ID']
    projects = df['Sample_Project']
    # send logging to the runlogger
    runlogger = setup_run_logger(run)
    ref_file, gtf = get_reference(run_dir)
    def sn_logger(sn_dict):
        if 'msg' in sn_dict:
            runlogger.info(sn_dict['msg'])
    output_dir = sample_sheet.get_output_dir()
    if output_dir:
        analysis =  output_dir
    else:
        analysis = run
    sm_config = {'run': run, 'samples': samples, 'projects': projects, 'ref': ref_file, 'gtf': gtf}
    opts = {
        'cores': 16,
        'nodes': 99,
        'local_cores': 4,
        'config': sm_config,
        'cluster_config': 'snakemake_cluster.json',
        'cluster': 'python slurm_submit.py',
        'cluster_status': 'python cluster_status.py',
        'printshellcmds': True,
        'printreason': True,
        #'cleanup_shadow': True,
        #'dryrun': True,
        'latency_wait': 60,
        #'touch': True,
        #'printdag': True,
        'log_handler': sn_logger
    }
    return opts

def run_snakemake(opts):
    ret = snakemake('Snakefile', **opts)
    return ret

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


def copy_log():
    # last line is the tail command
    try:
        logfile = logger.handlers[0].baseFilename
        lines = tail(logfile, 80000)
        # show max of 100 lines
        end = len(lines)
        max = 100
        if end > max:
            end = end - max
            lines = lines[-1:end:-1]
        else:
            lines = lines[::-1]
        with open(LOG_HTML, 'w') as f:
            f.write('<pre>')
            f.writelines(lines)
            f.write('</pre>')
    except Exception as e:
        logger.error('Error creating HTML log file: %s' % str(e))


def process_runs():
    '''
    Fills up the pool with runs with apply_async
    If anything is ready, result is placed in success_runs or failed_runs
    Then looks for more runs
    '''
    logger.info("Processing runs")

    run_dirs_tmp = find_runs(need_to_process)
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190506_A00794_0011_BHJTKKDMXX/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190506_A00794_0012_AHJTFFDMXX/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190506_NB502063_0322_AHGGVTBGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190506_NS500422_0818_AHGJN3BGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190508_NS500422_0820_AHGK5KBGXB/']
    run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/181011_NS500422_0732_AHMYF2BGX7/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190510_A00794_0015_BHJVWTDMXX/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190506_A00794_0012_AHJTFFDMXX/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/181011_NS500422_0732_AHMYF2BGX7/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190513_NB501677_0436_AHCVHGBGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190513_NB502063_0327_AHFL7YBGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190513_NB551608_0082_AHCFK3BGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190513_NS500422_0822_AHFFTMBGXB/']
    #run_dirs_tmp = ['/n/boslfs/INSTRUMENTS/illumina/190514_NB501677_0437_AHCGN5BGXB/']
    run_dirs = []
    for run in run_dirs_tmp:
        ss_path = get_sample_sheet_path(run)
        sample_sheet = SampleSheet(ss_path)
        sam_types = sample_sheet.get_sample_types()
        poly_A = sample_sheet.has_poly_A_index()
        if poly_A:
            continue
        for t, v in sam_types.items():
            if v == '10x single cell' or v == '10x single cell rna':
                run_dirs.append(run)
                break
    logger.info("Found %s runs: %s\n" % (len(run_dirs), json.dumps(run_dirs)))
    results = {}
    failed_runs = []
    success_runs = []
    for run_dir in run_dirs:
        #run_dir = run_dirs[0]

        logger.info("Current run  %s \n" % (json.dumps(run_dir)))
        run = os.path.basename(os.path.normpath(run_dir))
        opts = get_ody_snakemake_opts(run_dir)
        logger.info("Queueing odybcl2fastq cmd for %s:\n" % (run))
        success = run_snakemake(opts)

        if success:
            success_runs.append(run)
            status = 'success'
        else:
            failed_runs.append(run)
            status = 'failure'
            #failure_email(run, 'cmd', ret_code, 'std_out', 'std_err')

        logger.info(
            "Completed run: success %s and failures %s\n\n\n" %
            (json.dumps(success_runs), json.dumps(failed_runs))
        )


def main():
    try:
        logger.info("Starting ody10x processing")
        logger.info("Running with ")
        for k in ['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR', 'MOUNT_DIR', 'LOG_DIR', 'CONTROL_DIR']:
            logger.info("\t%s\t%s" % (k, config[k]))
        # run continuously
        #while True:
        # queue new runs for demultiplexing with bcl2fastq2
        process_runs()
        # check for any runs that started but never completed demultiplexing
        #notify_incomplete_runs()
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
