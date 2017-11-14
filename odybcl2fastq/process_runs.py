#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and run odybcl2fastq/parseargs.py on them

Created on  2017-11-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time
import logging
import subprocess
import json
from datetime import datetime
from multiprocessing import Pool
import odybcl2fastq.parseargs as parseargs
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage

SOURCE_DIR = '/n/boslfs/INSTRUMENTS/illumina/'
OUTPUT_DIR = '/n/regal/informatics/mportermahoney/odytest/'
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG_FILE = ROOT_DIR + '/odybcl2fastq.log'
PROCESSED_FILE = 'odybcl2fastq.processed'
COMPLETE_FILE = 'odybcl2fastq.complete'
INCOMPLETE_NOTIFIED_FILE = 'odybcl2fastq.incomplete_notified'
DAYS_TO_SEARCH = 7
INCOMPLETE_AFTER_DAYS = 0
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('Nov 1 2017', '%b %d %Y')
# TODO: seqprep looks for required file runParameters.xml or RunParameters.xml,
# do we need to check for that?
REQUIRED_FILES = ['SampleSheet.csv', 'InterOp/QMetricsOut.bin', 'InterOp/TileMetricsOut.bin', 'RunInfo.xml', 'RTAComplete.txt']
PROC_NUM = 3
FREQUENCY = 60

def setup_logging():
    # take level from env or INFO
    level = os.getenv('ODYBCL2FASTQ_LOGGING_LEVEL', logging.INFO)
    logging.basicConfig(
            filename= LOG_FILE,
            level=level,
            format='%(asctime)s %(message)s'
    )
    logging.getLogger().addHandler(logging.StreamHandler())

def failure_email(run, cmd, ret_code, std_out, std_err):
    log = parseargs.get_output_log(run)
    subject =  "Run Failed: %s" % run
    message = ("%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
            " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log))
    send_email(message, subject)

def send_email(message, subject):
    logging.warning(message)
    fromaddr = 'afreedman@fas.harvard.edu'
    toemaillist=['mportermahoney@g.harvard.edu']
    buildmessage(message, subject, None, fromaddr, toemaillist)

def touch(run_dir, file):
    # touch a processed file in the run_dir
    path = run_dir + file
    with open(path, 'w+'):
        os.utime(path, None)

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
    if ((now - m_time).days) < INCOMPLETE_AFTER_DAYS:
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
    dirs = glob.glob(SOURCE_DIR + '*/')
    runs = []
    for dir in dirs:
        if filter(dir):
            runs.append(dir)
    return runs

def runs_to_process():
    return find_runs(need_to_process)

def incomplete_runs():
    return find_runs(run_is_incomplete)

def get_odybcl2fastq_cmd(run_dir):
    run = os.path.basename(os.path.normpath(run_dir))
    params = {
        'runfolder': os.path.dirname(run_dir),
        'output-dir': OUTPUT_DIR + run,
        'sample-sheet': run_dir + 'SampleSheet.csv',
        'runinfoxml': run_dir + 'RunInfo.xml'
    }
    args = []
    opt_flag = '--'
    for opt, val in params.items():
        args.append(opt_flag +  opt)
        if val:
            args.append(val)
    return 'python ' + ROOT_DIR + '/odybcl2fastq/parseargs.py ' + ' '.join(args)

def run_odybcl2fastq(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    std_out, std_err = proc.communicate()
    return (proc.returncode, std_out, std_err, cmd)

def notify_incomplete_runs():
    run_dirs = incomplete_runs()
    run_dirs_str = "\n".join(run_dirs)
    message = "The following runs failed to complete %s or more days ago:\n\n%s" % (INCOMPLETE_AFTER_DAYS, run_dirs_str)
    send_email(message, 'Odybcl2fastq incomplete runs')
    for run in run_dirs:
        touch(run, INCOMPLETE_NOTIFIED_FILE)


def process_runs():
    runs_found = runs_to_process()
    proc_num = os.environ.get('ODYBCL2FASTQ_PROC_NUM', PROC_NUM)
    run_dirs = runs_found[:proc_num]
    logging.info("Found %s runs: %s\nprocessing first %s:\n%s\n" % (len(runs_found), json.dumps(runs_found), len(run_dirs),
        json.dumps(run_dirs)))
    pool = Pool(proc_num)
    results = {}
    for run_dir in run_dirs:
        run = os.path.basename(os.path.normpath(run_dir))
        cmd = get_odybcl2fastq_cmd(run_dir)
        logging.info("Queueing odybcl2fastq cmd for %s:\n%s\n" % (run, cmd))
        results[run] = pool.apply_async(run_odybcl2fastq, (cmd,))
        touch(run_dir, PROCESSED_FILE) # mark so it doesn't get reprocessed
    failed_runs = []
    success_runs = []
    for run, result in results.items():
        ret_code, std_out, std_err, cmd = result.get()
        logging.info("Odybcl2fastq for %s returned %i\n" % (run, ret_code))
        if ret_code == 0:
            success_runs.append(run)
        else:
            failed_runs.append(run)
            failure_email(run, cmd, ret_code, std_out, std_err)
        touch(run_dir, COMPLETE_FILE) # mark so it doesn't get reprocessed
    logging.info("Completed %i runs %i success %s and %i failures %s\n\n\n" %
            (len(results), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs)))

if __name__ == "__main__":
    try:
        setup_logging()
        # run continuously
        while True:
            process_runs()
            notify_incomplete_runs()
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logging.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
    except Exception as e:
        logging.exception(e)
