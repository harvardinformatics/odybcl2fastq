#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
parse arguments to run bcl2fastq, quality metrics analyses
and email reporting

Created on  2017-04-19

@author: Adam Freedman <adamfreedman@fas.harvard.edu>
@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import sys, os, stat
import logging
import json
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from collections import OrderedDict
import odybcl2fastq.util as util
from odybcl2fastq import config
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from subprocess import Popen, PIPE, STDOUT
from odybcl2fastq.status_db import StatusDB
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.run import Run
from odybcl2fastq.analysis import Analysis
from odybcl2fastq.qc.fastqc_runner import fastqc_runner

PROCESSED_FILE = 'odybcl2fastq.processed'
COMPLETE_FILE = 'odybcl2fastq.complete'
FINAL_DIR_PERMISSIONS = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH
FINAL_FILE_PERMISSIONS = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH
STORAGE_CAPACITY_WARN = 0.96
STORAGE_CAPACITY_ERROR = 0.99

logger = logging.getLogger('odybcl2fastq')


def initArgs():
    '''
    Setup arguments with parameterdef for bcl2fastq args
    other than the run and outpur dirs, check envs, parse
    commandline, return args
    '''

    parameterdefs = [
        {
            'name'      : 'TEST',
            'switches'  : ['-t', '--test'],
            'required'  : False,
            'help'      : 'run in test mode, log to std out, no not run cmd',
            'action'    : 'store_true',
        },
        {
            'name'      : 'NO_DEMULTIPLEX',
            'switches'  : ['--no-demultiplex'],
            'required'  : False,
            'help'      : 'run without demultiplexing, assumes fastq files are already in output dir',
            'action'    : 'store_true',
        },
        {
            'name'      : 'NO_POST_PROCESS',
            'switches'  : ['--no-post-process'],
            'required'  : False,
            'help'      : 'run without fastqc or other post processing',
            'action'    : 'store_true',
        },
        {
            'name'      : 'NO_FILE_COPY',
            'switches'  : ['--no-file-copy'],
            'required'  : False,
            'help'      : 'run without copying files to final',
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_MIN_LOG_LEVEL',
            'switches'  : ['--min-log-level'],
            'required'  : False,
            'help'      : 'logging level: none, fatal, error, warning, info, debug, trace',
            'default'   : 'INFO',
            'type'      : str,
        },
        {
            'name'      : 'BCL_PROC_THREADS',
            'switches'  : ['-p', '--processing-threads'],
            'required'  : False,
            'help'      : 'number threads to process demultiplexed data',
            'default'   : 8,
            'type'      : int,
        },
        {
            'name'      : 'BCL_ADAPTER_STRINGENCY',
            'switches'  : ['--adapter-stringency'],
            'required'  : False,
            'help'      : 'minimum match rate to trigger adapter masking/trimming',
            'default'   : 0.9,
            'type'      : float,
        },
        {
            'name'      : 'BCL_BARCODE_MISMATCHES',
            'switches'  : ['--barcode-mismatches'],
            'required'  : False,
            'help'      : 'number allowed mismatches per index',
            'default'   : 0,
            'type'      : int,
        },
        {
            'name'      : 'BCL_CREATE_INDEXREAD_FASTQ',
            'switches'  : ['--create-fastq-for-index-reads'],
            'required'  : False,
            'action'    : 'store_true',
            'help'      : 'writes fastq files for index read(s)',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_BCLS',
            'switches'  : ['--ignore-missing-bcls'],
            'required'  : False,
            'help'      : 'missing or corrupt bcl files are ignored',
            'default'   : True,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_FILTER',
            'switches'  : ['--ignore-missing-filter'],
            'required'  : False,
            'help'      : 'missing or corrupt filter files are ignored',
            'default'   : True,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_POSITIONS',
            'switches'  : ['--ignore-missing-positions'],
            'required'  : False,
            'help'      : 'missing or corrupt positions files are ignored',
            'default'   : True,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_MINIMUM_TRIMMED_READ_LENGTH',
            'switches'  : ['--minimum-trimmed-read-length'],
            'required'  : False,
            'help'      : 'minimum read length after adapter trimming',
            'default'   : 0,
            'type'      : int,
        },
        {
            'name'      : 'BCL_MASK_SHORT_ADAPTER_READS',
            'switches'  : ['--mask-short-adapter-reads'],
            'required'  : False,
            'help'      : 'controls adapter base masking if read below min length post trimming',
            'default'   : False,
            'type'      : int,
        },
        {
            'name'      : 'BCL_TILES',
            'switches'  : ['--tiles'],
            'required'  : False,
            'help'      : 'regex for tile selection',
            'default'   : False,
            'type'      : str,
        },
        {
            'name'      : 'BCL_WITH_FAILED_READS',
            'switches'  : ['--with-failed-reads'],
            'required'  : False,
            'help'      : 'include all clusters, including non-PF ones',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_WRITE_FASTQ_REVCOMP',
            'switches'  : ['--write-fastq-reverse-complement'],
            'required'  : False,
            'help'      : 'generate fastq files of reverse complements of actual data',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_NO_BGZF',
            'switches'  : ['--no-bgzf-compression'],
            'required'  : False,
            'help'      : 'turn off bgzf compression and use gzip instead for fastq files',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_NO_LANE_SPLITTING',
            'switches'  : ['--no-lane-splitting'],
            'required'  : False,
            'help'      : 'do not split fastq by lane',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_FIND_ADAPTERS_SLIDING_WINDOW',
            'switches'  : ['--find-adapters-with-sliding-window'],
            'required'  : False,
            'help'      : 'use simple sliding window to detec adapters, indels not handled',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_FASTQ_COMPRESSION_LEVEL',
            'switches'  : ['--fastq-compression-level'],
            'required'  : False,
            'help'      : 'use simple sliding window to detec adapters, indels not handled',
            'default'   : 4,
            'type'      : int,
            'choices'   : range(1, 10),
        },
        {
            'name'      : 'BCL_RUNFOLDER_DIR',
            'switches'  : ['-R', '--runfolder-dir'],
            'required'  : True,
            'help'      : 'path to run folder directory',
            'type'      : str,
        },
        {
            'name'      : 'BCL_OUTPUT_DIR',
            'switches'  : ['-o', '--output-dir'],
            'required'  : True,
            'help'      : 'path to demultiplexed output',
            'type'      : str,
        },
        {
            'name'      : 'BCL_SAMPLE_SHEET',
            'switches'  : ['--sample-sheet'],
            'required'  : True,
            'type'      : str,
            'help'      : 'path to sample sheet (if need to call custom sheets, e.g. for a within lane mixed index length run)',
        },
        {
            'name'      : 'RUNINFO_XML',
            'switches'  : '--runinfoxml',
            'required'  : False,
            'type'      : str,
            'help'      : 'path to runinfo xml file',
        },
        {
            'name'      : 'BCL_USE_BASES_MASK',
            'switches'  : '--use-bases-mask',
            'required'  : False,
            'default'   : False,
            'type'      : str,
            'help'      : 'base mask to use in demultiplexing, exp: y61,i8,i8,y14',
        },

    ]

    # Check for environment variable values
    # Set to 'default' if they are found
    for parameterdef in parameterdefs:
        if os.environ.get(parameterdef['name'], None) is not None:
            parameterdef['default'] = os.environ.get(parameterdef['name'])

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version='bcl2fastq2.2')

    # sets up dict with switches as keys and names as values
    switches_to_names = {}
    # Use the parameterdefs for the ArgumentParser
    for parameterdef in parameterdefs:
        switches = parameterdef.pop('switches')
        if not isinstance(switches, list):
            switches = [switches]

        # Gotta take it off for add_argument
        name = parameterdef.pop('name')
        parameterdef['dest'] = name
        if 'default' in parameterdef:
            parameterdef['help'] += '  [default: %s]' % parameterdef['default']
        parser.add_argument(*switches, **parameterdef)

        # Gotta put it back on for later
        parameterdef['name'] = name

        # add switch tuple as key, paramterdef name (= destination) as value
        if 'BCL' in parameterdef['name']:  # this allows non BCL things to be excluded from switches to names so don't get incorrectly added to cmd line arg
            switches_to_names[tuple(switches)] = parameterdef['name']
    args = parser.parse_args()
    return args, switches_to_names


def get_submissions(sample_sheet, instrument):
    subs = set()
    for key, row in sample_sheet['Data'].items():
        if row['Description']:
            subs.add(row['Description'])
    return list(subs)


def update_lims_db(run, sample_sheet, instrument):
    runlogger = logging.getLogger('run_logger')
    runlogger.info('Start db update for %s\n' % run)
    subs = get_submissions(sample_sheet, instrument)
    stdb = StatusDB()
    stdb.link_run_and_subs(run, subs)
    analysis = stdb.insert_analysis(run, ', '.join(subs))
    runlogger.info('End db update for %s\n' % analysis)


def copy_source_to_output(src_root, dest_root, sample_sheet, instrument):
    # copy important source files to the output dir they will then be moved to
    # final
    src_root += '/'
    dest_root += '/'
    if instrument == 'hiseq':
        run_params_file = 'runParameters.xml'
    else:  # nextseq
        run_params_file = 'RunParameters.xml'
    # get filename part of sample_sheet
    sample_sheet = sample_sheet.replace(src_root, '')
    files_to_copy = {
        'sample_sheet': sample_sheet,
        'run_info': 'RunInfo.xml',
        'run_params': run_params_file,
        'interop': 'InterOp'
    }
    for name, file in files_to_copy.items():
        new_file = file
        if name == 'sample_sheet':
            # needs to be called this for lims billing
            new_file = 'SampleSheet.csv'
        dest = dest_root  + new_file
        src = src_root + file
        util.copy(src, dest)


def copy_output_to_final(output_dir, run_folder, output_log):
    runlogger = logging.getLogger('run_logger')

    # determine dest_dir
    dest_dir = config.FINAL_DIR + run_folder
    # check size of output_dir
    cmd = 'du -s %s' % output_dir
    code, out, err = run_cmd(cmd)

    if code != 0:
        raise Exception('Could not check size of output_dir, files not copied to %s: %s' % (config.FINAL_DIR, cmd))
    output_space = int(out.split()[0])

    # check capacity of final dir
    cmd = 'df -P %s | grep  -v Filesystem' % config.FINAL_DIR
    code, out, err = run_cmd(cmd)
    if code != 0:
        raise Exception('Could not check capacity of %s, files not copied %s: %s' % (config.FINAL_DIR, output_dir, cmd))
    dest_space = out.split()
    tot_space = int(dest_space[1])
    used = int(dest_space[2])
    capacity = (used + output_space) / float(tot_space)
    storage_capacity_warn = float(os.getenv('ODY_STORAGE_CAPACITY_WARN', STORAGE_CAPACITY_WARN))
    storage_capacity_error = float(os.getenv('ODY_STORAGE_CAPACITY_ERROR', STORAGE_CAPACITY_ERROR))
    if capacity > storage_capacity_warn:
        message = '%s near capacity copying %s: %s, used: %s, tot: %s' % (config.FINAL_DIR, output_dir, output_space, used, tot_space)
        logger.warning(message)
        sent = buildmessage(message, 'NGSDATA is near capacity', [], config.EMAIL['from_email'], config.EMAIL['admin_email'])
    if capacity > storage_capacity_error:
        msg = 'Could not copy %s to  %s: %s' % (output_dir, config.FINAL_DIR, capacity)
        raise Exception(msg)
    runlogger.info('Copying %s: %s, to %s at capacity: %s' % (
        output_dir,
        output_space,
        dest_dir,
        capacity)
    )
    util.copy(output_dir, dest_dir)
    # change permissions on dest_dir
    util.chmod_rec(dest_dir, FINAL_DIR_PERMISSIONS, FINAL_FILE_PERMISSIONS)

def run_cmd(cmd):
    # run unix cmd, return out and error
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    out, err = proc.communicate()
    return (proc.returncode, out, err)

def analysis_runner(analysis, output_log, no_demultiplex=False):
    runlogger = logging.getLogger('run_logger')
    runlogger.info("***** START %s for %s *****\n\n" % (analysis.type, analysis.name))
    last_output = ''
    if no_demultiplex:
        message = 'run %s completed successfully\nsee logs here: %s\n' % (analysis.name, output_log)
        success = True
    else:
        runlogger.info('Launching %s...%s\n' % (analysis.type, analysis.cmd))
        code, last_output = analysis.execute(output_log)
        logger.info("***** END %s for %s *****\n\n" % (analysis.type, analysis.name))
        if code != 0:
            message = 'run %s failed\n see logs here: %s\n' % (analysis.name, output_log)
            success = False
        else:
            message = 'run %s completed successfully\nsee logs here: %s\n' % (analysis.name, output_log)
            success = True
    runlogger.info('message = %s' % message)
    return success, message + last_output

def process_runs(args=None, switches_to_names=None):
    ''' digest args and sample_sheet, call bcl2fastq or 10x demultiplexing
        args and switches are only passed in for the test
    '''
    if not args or not switches_to_names:
        args, switches_to_names = initArgs()
    logger.info('Processing runs from run folder %s, using sample sheet %s, to output dir %s' % (args.BCL_RUNFOLDER_DIR, args.BCL_SAMPLE_SHEET, args.BCL_OUTPUT_DIR))
    util.touch(args.BCL_RUNFOLDER_DIR + '/', PROCESSED_FILE)
    # create a run object with details about requested analysis
    run = Run(args, switches_to_names)
    runlogger = setup_run_logger(run.name, run.test)
    runlogger.info("***** START Odybcl2fastq *****\n\n")
    runlogger.info("Beginning to process run: %s\n args: %s\n" % (run.name, json.dumps(vars(args))))
    run.set_sample_sheet()

    # skip everything but billing if run folder flagged
    if os.path.exists(run.runfolder_dir + '/' + 'billing_only.txt'):
        runlogger.info("This run is flagged for billing only %s" % run.name)
        update_lims_db(run.name, run.sample_sheet.sections, run.instrument)
        return

    # extract indexing strategies from lanes, each strategy requires an analysis
    run.extract_basemasks()
    jobs_tot = len(run.mask_lists)
    if jobs_tot > 1:
        runlogger.info("This run contains different masks in the same lane and will require %i analysis jobs" % jobs_tot)
    job_cnt = 1
    # run an analysis per indexing strategy on run
    for mask, mask_list in run.mask_lists.items():
        output_suffix = ''
        # if more than one analysis cmd needed suffix output dir and sample sheet
        if jobs_tot > 1:
            output_suffix = mask.replace(',', '_')
        analysis = Analysis.create(run, mask, output_suffix)
        cmd = analysis.build_cmd(run.test)
        runlogger.info("\nJob %i of %i:" % (job_cnt, jobs_tot))
        if run.test:
            runlogger.info("Test run, command not run: %s" % cmd)
            success = True
            message = 'TEST'
        else:
            output_log = get_output_log(run.name)
            success, message = analysis_runner(analysis, output_log, run.no_demultiplex)
            summary_data = {}
            template = None
            # run folder will contain any suffix that was applied
            if success:
                if not run.no_post_process:
                    # update lims db
                    update_lims_db(analysis.name, analysis.sample_sheet.sections, run.instrument)
                    # run  qc, TODO: consider a seperate job for this
                    error_files, fastqc_err, fastqc_out = fastqc_runner(analysis.output_dir)
                    with open(output_log, 'a+') as f:
                        f.write('\n'.join(fastqc_out) + "\n\n")
                        f.write('\n'.join(fastqc_err) + "\n\n")
                # copy run files to final
                if not run.no_file_copy:
                    copy_source_to_output(
                        run.runfolder_dir,
                        analysis.output_dir,
                        analysis.sample_sheet_path,
                        run.instrument
                    )
                    run.fastq_checksum()
                    # copy output to final dest where users will access
                    copy_output_to_final(analysis.output_dir, analysis.name, output_log)
                # get data from analysis to put in the email
                summary_data = analysis.get_email_data()
                template = analysis.get_template()
                subject = 'Demultiplex Summary for ' + analysis.name
            else:
                subject = 'Run Failed: ' + analysis.name
            toemaillist = config.EMAIL['to_email']
            fromaddr = config.EMAIL['from_email']
            runlogger.info('Sending email summary to %s\n' % json.dumps(toemaillist))
            sent = buildmessage(message, subject, summary_data, template, fromaddr, toemaillist)
            runlogger.info('Email sent: %s\n' % str(sent))
            util.touch(run.output_dir + '/', COMPLETE_FILE)
        job_cnt += 1
    if success:
        ret_code = 0
        status = 'success'
        util.touch(run.runfolder_dir + '/', COMPLETE_FILE)
    else:
        # pass a special ret_code to avoid double email on error
        ret_code = 9
        status = 'failure'
    logger.info("odybcl2fastq for %s returned %s\n" % (run.name, status))
    runlogger.info("***** END Odybcl2fastq *****\n\n")
    return ret_code

def setup_run_logger(run_name, test):
    # take level from env or INFO
    runlogger = logging.getLogger('run_logger')
    level = logging.getLevelName(os.environ.get('ODYBCL2FASTQ_RUN_LOG_LEVEL', 'INFO'))
    runlogger.setLevel(level)
    handler = logging.FileHandler(get_output_log(run_name))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    handler.setLevel(level)
    handler.setFormatter(formatter)
    runlogger.addHandler(handler)
    return runlogger

def get_output_log(run):
    logdir = os.environ.get('ODYBCL2FASTQ_RUN_LOG_DIR', config.LOG_DIR)
    return os.path.join(logdir, run + '.log')

if __name__ == "__main__":
    try:
        sys.exit(process_runs())
    except Exception as e:
        logging.exception(e)
        raise
