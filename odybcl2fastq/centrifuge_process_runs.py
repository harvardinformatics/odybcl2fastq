#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
look for newly demultiplexed illumina runs and run centrifuge  on them

Created on  2018-07-05

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2018 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import os, glob, time, re
import logging
import subprocess
import json
from datetime import datetime
from multiprocessing import Pool
from odybcl2fastq import config
from odybcl2fastq import constants as const
import odybcl2fastq.util as util
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.run import COMPLETE_FILE as DEMULTIPLEX_COMPLETE_FILE
from odybcl2fastq.run import FINAL_DIR_PERMISSIONS, FINAL_FILE_PERMISSIONS

LOG_FILE = const.ROOT_DIR + 'centrifuge.log'
PROCESSED_FILE = 'centrifuge.processed'
COMPLETE_FILE = 'centrifuge.complete'
SKIP_FILE = 'centrifuge.skip'
FASTQLIST = 'centrifuge_fastqlist.txt'
DAYS_TO_SEARCH = 4
PROC_NUM = int(os.getenv('ODYBCL2FASTQ_PROC_NUM', 2))

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

def success_email(run, centrifuge_dir, cmd, ret_code, std_out, std_err):
    subject =  "Centrifuge Completed: %s" % run
    message = ("%s\nSee results: %s\ncmd: %s\nreturn code: %i" % (subject, centrifuge_dir, cmd, ret_code))
    send_email(message, subject)

def failure_email(run, cmd, ret_code, std_out, std_err):
    subject =  "Centrifuge Failed: %s" % run
    message = ("%s\ncmd: %s\nreturn code: %i\nstandard out: %s\nstandard"
            " error: %s\n" % (subject, cmd, ret_code, std_out, std_err))
    send_email(message, subject, 'centrifuge_admin_email')

def send_email(message, subject, to_email = None):
    logging.warning(message)
    fromaddr = config.EMAIL['centrifuge_from_email']
    if to_email:
        toemaillist = config.EMAIL[to_email]
    else:
        toemaillist=config.EMAIL['centrifuge_to_email']
    buildmessage(message, subject, None, fromaddr, toemaillist)

def need_to_process(dir):
    run = dir.split('/')[-2]
    now = datetime.now()
    m_time = datetime.fromtimestamp(os.stat(dir).st_mtime)
    # filter out if modified outside or search window
    if ((now - m_time).days) > DAYS_TO_SEARCH:
        return False
    # filter out if tagged as processed
    if os.path.isfile(dir + PROCESSED_FILE):
        return False
    # filter out if tagged as skip
    if os.path.isfile(dir + SKIP_FILE):
        return False
    # filter out if run never completed to get transfered to ngsdata
    if not os.path.exists(config.FINAL_DIR + run):
        return False
    # filter out if run never completed demultiplexing
    if not os.path.exists(dir + DEMULTIPLEX_COMPLETE_FILE ):
        return False
    return True

def find_runs(filter):
    # get all subdirectories
    dirs = sorted(glob.glob(config.OUTPUT_DIR + '*/'))
    runs = []
    for dir in dirs:
        if filter(dir):
            runs.append(dir)
    return runs

def get_fastq_files(dir):
    # check if we are limiting to a list of fastq files
    if os.path.isfile(dir + FASTQLIST):
        file_lst = []
        with open(dir + FASTQLIST) as f:
            file_lst = [line.strip() for line in f if os.path.isfile(line.strip())]
        if file_lst:
            return file_lst
    # get fastq files except undetermined
    sample_proj_path = '%s/*/*.fastq.gz' % dir
    file_path = '%s/*.fastq.gz' % dir
    file_lst = glob.glob(sample_proj_path)
    file_lst.extend(glob.glob(file_path))
    # exclude undetermined files
    file_lst = [path for path in file_lst if not re.match(r'.*/Undetermined.*\.fastq\.gz', path)]
    return file_lst

def group_fastq_files(run_dir, files):
    # fastq files must be grouped by sample_name/sample then read then lane
    grps = {}
    for path in files:
        path_lst = path
        path_lst = path_lst.replace(run_dir, '').split('/')
        grp = ''
        if len(path_lst) > 1:
            grp = path_lst.pop(0) + '_'
        if len(path_lst) > 1:
            raise Exception('fastq path deeper than expected: %s' % path)
        path_lst = path_lst[0].replace('.fastq.gz', '').split('_')

        path_lst.pop() # remove the 001
        read = int(path_lst.pop().replace('R', ''))
        lane = path_lst.pop()
        if 'L' in lane:
            lane = int(lane.replace('L', ''))
        else: # nextseq does not have lanes so group all in 1
            lane = 1
        path_lst.pop() # remove the S
        sample = '_'.join(path_lst)
        grp += sample
        if grp not in grps:
            grps[grp] = {}
        if read not in grps[grp]:
            grps[grp][read] = {}
        grps[grp][read][lane] = path
    return grps

def get_centrifuge_cmd(run_dir, grp, read1, read2):
    run = os.path.basename(os.path.normpath(run_dir))
    # create centrifuge dir
    centrifuge_dir = run_dir + 'centrifuge'
    outfile = centrifuge_dir + '/' + grp + '.html'
    if not os.path.exists(centrifuge_dir):
        subprocess.call('mkdir %s' % (centrifuge_dir) ,shell=True)
    script = config.CENTRIFUGE_DIR + 'centrifuge.sh'
    cmd_lst = ['bash -e', script, ','.join(read1)]
    if read2:
        cmd_lst.append(','.join(read2))
    else:
        cmd_lst.append('None')
    cmd_lst.append(outfile)
    cmd =  ' '.join(cmd_lst)
    return cmd, outfile

def run_centrifuge(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    std_out, std_err = proc.communicate()
    return (proc.returncode, std_out, std_err, cmd)

def process_runs(pool, proc_num):
    runs_found = find_runs(need_to_process)
    if runs_found:
        run_dir = runs_found[0]
        util.touch(run_dir, PROCESSED_FILE)
        logging.info("Found %s runs: %s\nprocessing the first:\n%s\n" % (len(runs_found), json.dumps(runs_found),
            run_dir))
        results = {}
        run = os.path.basename(os.path.normpath(run_dir))
        logging.info("Centrifuge processing %s\n" % (run))
        # we run centrifuge per sample
        fastq_files = get_fastq_files(run_dir)
        grps = group_fastq_files(run_dir, fastq_files)
        logging.info("Found %s samples to run on centrifuge\n" % (len(grps)))
        for sample, files in grps.items():
            read1 = [files[1][i] for i in sorted(files[1])]
            if 2 in files:
                read2 = [files[2][i] for i in sorted(files[2])]
            else:
                read2 = None
            cmd, output_dir = get_centrifuge_cmd(run_dir, sample, read1, read2)
            centrifuge_dir = run_dir + 'centrifuge'
            logging.info("Queueing centrifuge cmd for %s:\n%s\n" % (sample, cmd))
            results[sample] = pool.apply_async(run_centrifuge, (cmd,))
        failed_samples = []
        success_samples = []
        for sample, result in results.items():
            ret_code, std_out, std_err, cmd = result.get()
            output = std_out + std_err
            if ret_code == 0:
                success_samples.append(sample)
                message = 'sample %s completed successfully\nsee logs here: %s\n' % (sample, output)
                status = 'success'
            else:
                failed_samples.append(sample)
                message = 'sample %s failed\n see logs here: %s\n' % (sample, output)
                status = 'failure'
                # failures from bcl2fastq will be emailed from inner job
                # inner job passes ret code 9 on fail from bcl2fastq
                # only email from outer job if error is from inner job itself
                # not the bcl2fastq subprocess
            logging.info('message = %s' % message)

        if failed_samples:
            failure_email(run, cmd, ret_code, std_out, std_err)
            # success or failure of individual run will be logged from run.py to capture
            # manual runs for the status log
        else:
            # copy html files to final dir
            dest_dir = config.FINAL_DIR + run + '/centrifuge/'
            centrifuge_dir = run_dir + 'centrifuge'
            html_files = glob.glob(centrifuge_dir + '/*.html')
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            for hf in html_files:
                html_name = hf.split('/')[-1]
                util.copy(hf, dest_dir + html_name)
            util.chmod_rec(dest_dir, FINAL_DIR_PERMISSIONS, FINAL_FILE_PERMISSIONS)
            util.touch(run_dir, COMPLETE_FILE)
            success_email(run, centrifuge_dir, cmd, ret_code, std_out, std_err)
        logging.info("Completed centrifuge for run %s with %i samples %i success %s and %i failures %s\n\n\n" %
                (run, len(results), len(success_samples), json.dumps(success_samples), len(failed_samples), json.dumps(failed_samples)))

if __name__ == "__main__":
    try:
        setup_logging()
        proc_num = PROC_NUM
        # create pool and call process_runs to apply_async jobs
        pool = Pool(proc_num)
        # run continuously
        while True:
            # queue new runs for demultiplexing with bcl2fastq2
            process_runs(pool, proc_num)
            # wait before checking for more runs to process
            frequency = os.getenv('ODYBCL2FASTQ_FREQUENCY', FREQUENCY)
            if frequency != FREQUENCY:
                logging.info("Frequency is not default: %i\n" % frequency)
            time.sleep(frequency)
        pool.close()
    except Exception as e:
        logging.exception(e)
        send_email(str(e), 'Odybcl2fastq centrifuge exception', 'admin_email')
