#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and enter then into bauer db

Created on  2017-11-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time, sys
import json
from datetime import datetime
from odybcl2fastq import config, initLogger
import odybcl2fastq.util as util
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.bauer_db import BauerDB
from odybcl2fastq.run import get_output_log

PROCESSED_FILE = 'bauer.processed'
COMPLETE_FILE = 'bauer.complete'
INCOMPLETE_NOTIFIED_FILE = 'bauer.incomplete_notified'
INCOMPLETE_AFTER_DAYS = 1
DAYS_TO_SEARCH = 7
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('May 10 2018', '%b %d %Y')
REQUIRED_FILES = ['SampleSheet.csv', 'RunInfo.xml']
PROC_NUM = 1

# Setup logger
logger = initLogger('load_runs')

# wait before checking for more runs to process
DEFAULT_SLEEP = 60
SLEEPYTIME = int(os.getenv('ODYBCL2FASTQ_SLEEP', DEFAULT_SLEEP))
if SLEEPYTIME != DEFAULT_SLEEP:
    logger.info("Frequency is not default: %i\n" % SLEEPYTIME)


def failure_email(run, cmd, ret_code, std_out, std_err):
    log = get_output_log(run)
    subject = "Run DB insert Failed: %s" % run
    message = (
        "%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
        " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log)
    )
    send_email(message, subject)


def send_email(message, subject):
    logger.warning(message)
    fromaddr = config.EMAIL['from_email']
    toemaillist = config.EMAIL['to_email']
    buildmessage(message, subject, None, fromaddr, toemaillist)


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
    if os.path.isfile(dir + '/' + COMPLETE_FILE):
        return False
    # filter out if never tagged for processing
    if not os.path.isfile(dir + '/' + PROCESSED_FILE):
        return False
    # filter out already notified
    if os.path.isfile(dir + '/' + INCOMPLETE_NOTIFIED_FILE):
        return False
    return True


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
    if os.path.isfile(dir + '/' + PROCESSED_FILE):
        return False
    # filter out if any required files are missing
    for req in REQUIRED_FILES:
        if not os.path.exists(dir + '/' + req):
            return False
    return True


def find_runs(filter):
    # get all subdirectories
    searchpath = os.path.join(config.SOURCE_DIR, '*')
    dirs = sorted(glob.glob(searchpath))
    runs = []
    for d in dirs:
        if filter(d):
            runs.append(d)
    logger.info('Found %d %s runs from %s' % (len(runs), filter.__name__, searchpath))
    return runs


def get_sample_sheet_path(run_dir):
    # set default
    sample_sheet_path = run_dir + '/SampleSheet.csv'
    return sample_sheet_path


def notify_incomplete_runs():
    run_dirs = find_runs(run_is_incomplete)
    run_dirs_str = "\n".join(run_dirs)
    if run_dirs:
        message = "The following runs failed be entered into bauer db %s or more days ago:\n\n%s" % (INCOMPLETE_AFTER_DAYS, run_dirs_str)
        send_email(message, 'BauerDB incomplete runs')
        for run in run_dirs:
            util.touch(run + '/', INCOMPLETE_NOTIFIED_FILE)


def load_runs(proc_num):
    runs_found = find_runs(need_to_process)
    run_dirs = runs_found[:proc_num]
    if run_dirs:
        logger.info(
            "Found %s runs: %s\nprocessing first %s:\n%s\n" % (len(runs_found), json.dumps(runs_found), len(run_dirs), json.dumps(run_dirs)))
        results = {}
        for run_dir in run_dirs:
            util.touch(run_dir + '/', PROCESSED_FILE)
            run = os.path.basename(os.path.normpath(run_dir))
            sample_sheet_path = get_sample_sheet_path(run_dir)
            bauer = BauerDB(sample_sheet_path)
            logger.info("Loading run into bauer db:\n%s\n" % (run))
            results[run] = bauer.insert_run(run_dir + '/')
        failed_runs = []
        success_runs = []
        for run, result in results.items():
            if result:
                success_runs.append(run)
                util.touch(run_dir + '/', COMPLETE_FILE)
            else:
                failed_runs.append(run)

            # success or failure of individual run will be logged from run.py to capture
            # manual runs for the status log
        logger.info(
            "Completed %i runs %i success %s and %i failures %s\n\n\n" %
            (len(results), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs))
        )


def main():
    try:
        proc_num = os.getenv('ODYBCL2FASTQ_PROC_NUM', PROC_NUM)
        nodaemon = len(sys.argv) > 1 and sys.argv[1] == '--no-daemon'
        if nodaemon:
            logger.info('--no-daemon set; running a single pass')
        # run continuously
        while True:
            # search for new runs
            load_runs(proc_num)
            notify_incomplete_runs()

            # If --no-daemon argument is present, break
            if nodaemon:
                break

            time.sleep(SLEEPYTIME)

    except Exception as e:
        logger.exception(e)
        send_email(str(e), 'Odybcl2fastq load_runs.py exception')


if __name__ == "__main__":
    sys.exit(main())
