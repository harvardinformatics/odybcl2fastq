#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and enter then into bauer db

Created on  2017-11-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time
import logging
import json
from datetime import datetime
from odybcl2fastq import config
from odybcl2fastq import constants as const
import odybcl2fastq.util as util
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.bauer_db import BauerDB

LOG_FILE = const.ROOT_DIR + 'db.log'
PROCESSED_FILE = 'bauer.processed'
COMPLETE_FILE = 'bauer.complete'
DAYS_TO_SEARCH = 7
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('May 07 2018', '%b %d %Y')
REQUIRED_FILES = ['SampleSheet.csv', 'RunInfo.xml']
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
    subject =  "Run DB insert Failed: %s" % run
    message = ("%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
            " error: %s\nsee log: %s\n" % (subject, cmd, ret_code, std_out, std_err, log))
    send_email(message, subject)

def send_email(message, subject):
    logging.warning(message)
    fromaddr = config.EMAIL['from_email']
    toemaillist=config.EMAIL['to_email']
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
    # filter out if any required files are missing
    for req in REQUIRED_FILES:
        if not os.path.exists(dir + req):
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

def load_runs(proc_num):
    runs_found = find_runs(need_to_process)
    run_dirs = runs_found[:proc_num]
    if run_dirs:
        logging.info("Found %s runs: %s\nprocessing first %s:\n%s\n" % (len(runs_found), json.dumps(runs_found), len(run_dirs),
            json.dumps(run_dirs)))
        results = {}
        for run_dir in run_dirs:
            run = os.path.basename(os.path.normpath(run_dir))
            sample_sheet_path = get_sample_sheet_path(run_dir)
            bauer = BauerDB(sample_sheet_path)
            logging.info("Loading run into bauer db:\n%s\n" % (run))
            results[run] = bauer.insert_run(run_dir)
        failed_runs = []
        success_runs = []
        for run, result in results.items():
            if result:
                success_runs.append(run)
                status = 'success'
            else:
                failed_runs.append(run)
                status = 'failure'

            # success or failure of individual run will be logged from run.py to capture
            # manual runs for the status log
        logging.info("Completed %i runs %i success %s and %i failures %s\n\n\n" %
                (len(results), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs)))

if __name__ == "__main__":
    try:
        proc_num = os.getenv('ODYBCL2FASTQ_PROC_NUM', PROC_NUM)
        setup_logging()
        # run continuously
        while True:
            # search for new runs
            load_runs(proc_num)
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logging.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
    except Exception as e:
        logging.exception(e)
