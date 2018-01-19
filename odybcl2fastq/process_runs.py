#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and run odybcl2fastq/odybcl2fastq.py on them

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
from odybcl2fastq import config
from odybcl2fastq import constants as const
from odybcl2fastq import run as ody_run
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage

LOG_FILE = const.ROOT_DIR + 'odybcl2fastq.log'
LOG_HTML = config.FINAL_DIR + 'odybcl2fastq_log.html'
PROCESSED_FILE = 'odybcl2fastq.processed'
COMPLETE_FILE = 'odybcl2fastq.complete'
INCOMPLETE_NOTIFIED_FILE = 'odybcl2fastq.incomplete_notified'
DAYS_TO_SEARCH = 7
INCOMPLETE_AFTER_DAYS = 2
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('Jan 15 2017', '%b %d %Y')
REQUIRED_FILES = ['InterOp/QMetricsOut.bin', 'InterOp/TileMetricsOut.bin', 'RunInfo.xml', 'RTAComplete.txt']
PROC_NUM = 2
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
    log = ody_run.get_output_log(run)
    subject =  "Run Failed: %s" % run
    message = ("%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
            " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log))
    send_email(message, subject)

def send_email(message, subject):
    logging.warning(message)
    fromaddr = config.EMAIL['from_email']
    toemaillist=config.EMAIL['admin_email']
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

def get_odybcl2fastq_cmd(run_dir):
    run = os.path.basename(os.path.normpath(run_dir))
    params = {
        'runfolder': os.path.dirname(run_dir),
        'output-dir': config.OUTPUT_DIR + run,
        'sample-sheet': run_dir + 'SampleSheet.csv',
        'runinfoxml': run_dir + 'RunInfo.xml'
    }
    args = []
    opt_flag = '--'
    for opt, val in params.items():
        args.append(opt_flag +  opt)
        if val:
            args.append(val)
    return 'python ' + const.APP_DIR + 'run.py ' + ' '.join(args)

def run_odybcl2fastq(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    std_out, std_err = proc.communicate()
    return (proc.returncode, std_out, std_err, cmd)

def notify_incomplete_runs():
    run_dirs = find_runs(run_is_incomplete)
    run_dirs_str = "\n".join(run_dirs)
    if run_dirs:
        message = "The following runs failed to complete %s or more days ago:\n\n%s" % (INCOMPLETE_AFTER_DAYS, run_dirs_str)
        send_email(message, 'Odybcl2fastq incomplete runs')
        for run in run_dirs:
            touch(run, INCOMPLETE_NOTIFIED_FILE)
def tail(f, n):
    proc = subprocess.Popen("tail -n %i %s | grep returned" % (n, f), shell=True, stdout=subprocess.PIPE)
    std_out, std_err = proc.communicate()
    return std_out.splitlines(True)

def copy_log():
    # last line is the tail command
    lines = tail(LOG_FILE, 5000)[:-1]
    # show max of 100 lines
    end = len(lines)
    if end > 100:
        end = 100
    lines = lines[end::-1]
    with open(LOG_HTML, 'w') as f:
        f.write('<pre>')
        f.writelines(lines)
        f.write('</pre>')

def process_runs(pool, proc_num):
    runs_found = find_runs(need_to_process)
    run_dirs = runs_found[:proc_num]
    if run_dirs:
        logging.info("Found %s runs: %s\nprocessing first %s:\n%s\n" % (len(runs_found), json.dumps(runs_found), len(run_dirs),
            json.dumps(run_dirs)))
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
            if ret_code == 0:
                success_runs.append(run)
                touch(run_dir, COMPLETE_FILE)
                status = 'success'
            else:
                failed_runs.append(run)
                status = 'failure'
            logging.info("Odybcl2fastq for %s returned %s\n" % (run, status))
        if failed_runs:
            send_email(('The following runs failed: ' + json.dumps(failed_runs)), 'Odybcl2fastq failed runs')
        logging.info("Completed %i runs %i success %s and %i failures %s\n\n\n" %
                (len(results), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs)))
    copy_log()

if __name__ == "__main__":
    try:
        setup_logging()
        proc_num = os.getenv('ODYBCL2FASTQ_PROC_NUM', PROC_NUM)
        # create pool and call process_runs to apply_async jobs
        pool = Pool(proc_num)
        # run continuously
        while True:
            # queue new runs for demultiplexing with bcl2fastq2
            process_runs(pool, proc_num)
            # check for any runs that started but never completed demultiplexing
            notify_incomplete_runs()
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logging.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
        pool.close()
    except Exception as e:
        logging.exception(e)
