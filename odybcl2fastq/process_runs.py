#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly completed illumina runs and run odybcl2fastq/odybcl2fastq.py on them

Created on  2017-11-01

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
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

LOG_HTML = config.FINAL_DIR + 'odybcl2fastq_log.html'
PROCESSED_FILE = ody_run.PROCESSED_FILE
COMPLETE_FILE = ody_run.COMPLETE_FILE
SKIP_FILE = 'odybcl2fastq.skip'
INCOMPLETE_NOTIFIED_FILE = 'odybcl2fastq.incomplete_notified'
DAYS_TO_SEARCH = 7
INCOMPLETE_AFTER_DAYS = 2
# a hardcoded date not to search before
# this will be helpful in transitioning from seqprep to odybcl2fastq
SEARCH_AFTER_DATE = datetime.strptime('Jan 15 2017', '%b %d %Y')
REQUIRED_FILES = ['InterOp/QMetricsOut.bin', 'InterOp/TileMetricsOut.bin', 'RunInfo.xml', 'RTAComplete.txt']
PROC_NUM = int(os.getenv('ODYBCL2FASTQ_PROC_NUM', 2))

FREQUENCY = 60

logger = logging.getLogger('odybcl2fastq')


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


def get_odybcl2fastq_cmd(run_dir):
    run = os.path.basename(os.path.normpath(run_dir))
    params = {
        'runfolder': os.path.dirname(run_dir),
        'output-dir': config.OUTPUT_DIR + run,
        'sample-sheet': get_sample_sheet_path(run_dir),
        'runinfoxml': run_dir + 'RunInfo.xml'
    }
    args = []
    opt_flag = '--'
    for opt, val in params.items():
        args.append(opt_flag + opt)
        if val:
            args.append(val)
    return 'python ' + const.APP_DIR + 'run.py ' + ' '.join(args)


def run_odybcl2fastq(cmd):
    proc = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
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


def process_runs(pool, proc_num):
    '''
    Fills up the pool with runs with apply_async
    If anything is ready, result is placed in success_runs or failed_runs
    Then looks for more runs
    '''
    logger.info("Processing runs")

    # Filtering out 10x runs to prevent the pipeline from processing
    # TODO: remove this when we can move this pipeline to snakemake
    run_dirs_all = find_runs(need_to_process)
    logger.info("Found %s runs: %s\n" % (len(run_dirs_all), json.dumps(run_dirs_all)))
    run_dirs = []
    for run in run_dirs_all:
        ss_path = get_sample_sheet_path(run)
        sample_sheet = SampleSheet(ss_path)
        sam_types = sample_sheet.get_sample_types()
        types10x = [
            '10x single cell',
            '10x genome',
            '10x single cell rna',
            '10x singel cell atac'
        ]
        is_10x = False
        for t, v in sam_types.items():
            if v in types10x:
                is_10x = True
                break
        if not is_10x:
            run_dirs.append(run)
    logger.info("Found %s runs that are not 10x: %s\n" % (len(run_dirs), json.dumps(run_dirs)))
    results = {}
    failed_runs = []
    success_runs = []
    queued_runs = {}

    while len(run_dirs) > 0 or len(results) > 0:
        logger.info("Current runs found: %s Runs pending results: %s Runs queued: %s\n" % (json.dumps(run_dirs), json.dumps(results.keys()), json.dumps(queued_runs.keys())))
        while len(run_dirs) > 0:
            run_dir = run_dirs.pop()
            run = os.path.basename(os.path.normpath(run_dir))
            cmd = get_odybcl2fastq_cmd(run_dir)
            logger.info("Queueing odybcl2fastq cmd for %s:\n%s\n" % (run, cmd))
            results[run] = pool.apply_async(run_odybcl2fastq, (cmd,))
            queued_runs[run] = 1

        for run, result in results.items():
            if result.ready():
                ret_code, std_out, std_err, cmd = result.get()
                if ret_code == 0:
                    success_runs.append(run)
                    status = 'success'
                else:
                    failed_runs.append(run)
                    status = 'failure'
                    # failures from bcl2fastq will be emailed from inner job
                    # inner job passes ret code 9 on fail from bcl2fastq
                    # only email from outer job if error is from inner job itself
                    # not the bcl2fastq subprocess
                    if ret_code != 9:
                        failure_email(run, cmd, ret_code, std_out, std_err)

                # Clear this result from the dict
                del results[run]
                del queued_runs[run]
        logger.info("After checking results, runs found: %s Runs pending results: %s Runs queued %s\n" % (json.dumps(run_dirs), json.dumps(results.keys()), json.dumps(queued_runs.keys())))

        sleep(10)
        new_run_dirs = find_runs(need_to_process)
        for new_run_dir in new_run_dirs:
            new_run = os.path.basename(os.path.normpath(new_run_dir))
            if new_run not in queued_runs:
                run_dirs.append(new_run_dir)

        if len(run_dirs) > 0:
            logger.info("Found %s more runs: %s\n" % (len(run_dirs), json.dumps(run_dirs)))


    logger.info(
        "Completed %i runs %i success %s and %i failures %s\n\n\n" %
        (len(results), len(success_runs), json.dumps(success_runs), len(failed_runs), json.dumps(failed_runs))
    )
    copy_log()


def main():
    try:
        proc_num = PROC_NUM
        # create pool and call process_runs to apply_async jobs
        pool = Pool(proc_num)
        logger.info("Starting odybcl2fastq processing")
        logger.info("Running with ")
        for k in ['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR', 'MOUNT_DIR', 'LOG_DIR', 'CONTROL_DIR']:
            logger.info("\t%s\t%s" % (k, config[k]))
        # run continuously
        while True:
            # queue new runs for demultiplexing with bcl2fastq2
            process_runs(pool, proc_num)
            # check for any runs that started but never completed demultiplexing
            notify_incomplete_runs()
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logger.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
        pool.close()
    except Exception as e:
        logging.exception(e)
        send_email(str(e), 'Odybcl2fastq exception')
        return 1


if __name__ == "__main__":
    sys.exit(main())
