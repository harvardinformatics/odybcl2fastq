#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
parse arguments to run bcl2fastq, quality metrics analyses
and email reporting

Created on  2017-04-19

@author: Adam Freedman <adamfreedman@fas.harvard.edu>
@copyright: 2017 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import sys, os, traceback, stat
import logging
import json
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import odybcl2fastq.util as util
from odybcl2fastq import config
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.parsers import parse_stats
from subprocess import Popen,PIPE
from odybcl2fastq.status_db import StatusDB
import odybcl2fastq.parsers.parse_sample_sheet as ss
from odybcl2fastq.qc.fastqc_runner import fastqc_runner
from tests.compare_fastq import compare_fastq

FINAL_DIR_PERMISSIONS = stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IXGRP|stat.S_IROTH|stat.S_IXOTH
FINAL_FILE_PERMISSIONS = stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH
INDROP_FILE = 'indrop.txt'

def initArgs():
    '''
    Setup arguments with parameterdef for bcl2fastq args
    other than the run and outpur dirs, check envs, parse
    commandline, return args
    '''

    parameterdefs = [
        {
            'name'      : 'TEST',
            'switches'  : ['-t','--test'],
            'required'  : False,
            'help'      : 'run in test mode, log to std out, no not run cmd',
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_PROC_THREADS',
            'switches'  : ['-p','--processing-threads'],
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
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_FILTER',
            'switches'  : ['--ignore-missing-filter'],
            'required'  : False,
            'help'      : 'missing or corrupt filter files are ignored',
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_POSITIONS',
            'switches'  : ['--ignore-missing-positions'],
            'required'  : False,
            'help'      : 'missing or corrupt positions files are ignored',
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
            'default'   : 22,
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
            'choices'   : range(1,10),
        },
        {
            'name'      : 'BCL_RUNFOLDER_DIR',
            'switches'  : ['-R','--runfolder-dir'],
            'required'  : True,
            'help'      : 'path to run folder directory',
            'type'      : str,
        },
        {
            'name'      : 'BCL_OUTPUT_DIR',
            'switches'  : ['-o','--output-dir'],
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
        if os.environ.get(parameterdef['name'],None) is not None:
            parameterdef['default'] = os.environ.get(parameterdef['name'])

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version='bcl2fastq2.19')

    # sets up dict with switches as keys and names as values
    switches_to_names={}
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
        parser.add_argument(*switches,**parameterdef)

        # Gotta put it back on for later
        parameterdef['name'] = name

        # add switch tuple as key, paramterdef name (= destination) as value
        if 'BCL' in parameterdef['name']: # this allows non BCL things to be excluded from switches to names so don't get incorrectly added to cmd line arg
            switches_to_names[tuple(switches)]=parameterdef['name']
    args = parser.parse_args()
    return args, switches_to_names

def get_submissions(sample_sheet, instrument):
    subs = set()
    for key, row in sample_sheet['Data'].items():
        if row['Description']:
            subs.add(row['Description'])
    return list(subs)

def update_lims_db(run, sample_sheet, instrument):
    logging.info('Start db update for %s\n' % run)
    subs = get_submissions(sample_sheet, instrument)
    stdb = StatusDB()
    analysis = stdb.insert_analysis(run, ', '.join(subs))
    stdb.link_run_and_subs(run, subs)
    logging.info('End db update for %s\n' % analysis)

def copy_source_to_output(src_root, dest_root, sample_sheet, instrument):
    # copy important source files to the output dir they will then be moved to
    # final
    src_root += '/'
    dest_root += '/'
    if instrument == 'hiseq':
        run_params_file = 'runParameters.xml'
    else: # nextseq
        run_params_file = 'RunParameters.xml'
    # get filename part of sample_sheet
    sample_sheet = sample_sheet.replace(src_root, '')
    files_to_copy = [
            sample_sheet,
            'RunInfo.xml',
            run_params_file,
            'InterOp'
    ]
    for file in files_to_copy:
        dest = dest_root  + file
        src = src_root + file
        util.copy(src, dest)

def copy_output_to_final(output_dir, run_folder, suffix):
    # determine dest_dir
    dest_dir = config.FINAL_DIR + run_folder
    if suffix: # runs with multiple indexing strategies have a subdir
        dest_dir += '/' + suffix
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
    if capacity > 0.8:
        logging.warning('%s near capacity copying %s: %s, used: %s, tot: %s' % (config.FINAL_DIR, output_dir, output_space, used, tot_space))
    if capacity > 0.9:
        msg = 'Could not copy %s to  %s: %s' % (output_dir, config.FINAL_DIR, capacity)
        raise Exception(msg)
    logging.info('Copying %s: %s, to %s at capacity: %s' % (output_dir, output_space,
        dest_dir, capacity))
    util.copy(output_dir, dest_dir)
    # change permissions on dest_dir
    util.chmod_rec(dest_dir, FINAL_DIR_PERMISSIONS, FINAL_FILE_PERMISSIONS)

def run_cmd(cmd):
    # run unix cmd, return out and error
    proc = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    out, err = proc.communicate()
    return (proc.returncode, out, err)

def bcl2fastq_build_cmd(args, switches_to_names, mask_list, instrument, run_type):
    argdict = vars(args)
    mask_switch = '--use-bases-mask'
    # each mask should be prefaced by the switch
    mask_opt = mask_switch + ' ' + (' ' + mask_switch + ' ').join(mask_list)
    cmdstrings=['bcl2fastq', mask_opt]
    # keeps consistent order of writing
    switch_list=switches_to_names.keys()
    switch_list.sort()
    for switches in switch_list:
        switch=[switch for switch in switches if '--' in switch][0]
        argvalue=str(argdict[switches_to_names[switches]])
        # the bit below prevents boolean flags from having values in the cmd
        if argvalue != 'False':
            if argvalue == 'True':
                cmdstrings.append(switch)
            else:
                cmdstrings.append(' '.join([switch,argvalue]))
    if instrument == 'nextseq':
        cmdstrings.append('--no-lane-splitting')
    if run_type == 'indrop':
        cmdstrings.append('--create-fastq-for-index-reads')
    cmdstring=' '.join(cmdstrings)
    return cmdstring

def bcl2fastq_runner(cmd,args):
    logging.info("***** START bcl2fastq *****\n\n")
    run = os.path.basename(args.BCL_RUNFOLDER_DIR)
    output_log = get_output_log(run)
    code, demult_out, demult_err = run_cmd(cmd)
    # append to output to log for the run
    with open(output_log, 'a+') as f:
        f.write(demult_err + "\n\n")
    logging.info("***** END bcl2fastq *****\n\n")
    if code!=0:
        message = 'run %s failed\n see logs here: %s\n%s\n' % (run, output_log,
                demult_err)
        success = False
    else:
        message = 'run %s completed successfully\nsee logs here: %s\n' % (run, output_log)
        success = True
    return success, message

def write_new_sample_sheet(new_samples, sample_sheet, output_suffix):
    new_sample_sheet = sample_sheet.replace('.csv', ('_' + output_suffix + '.csv'))
    input = open(sample_sheet, 'r')
    output = open(new_sample_sheet, 'wb')
    for line in input:
        if not line.startswith('[Data]'):
            output.write(line)
        else:
            output.write(line)
            break
    # print data headers
    output.write(next(input))
    # write new samples to sheet
    new_lines = [(','.join(row.values()) + "\r\n") for row in new_samples]
    output.writelines(new_lines)
    output.close()
    input.close()
    return new_sample_sheet

def get_run_type(run_dir):
    # flag files will be used for some special run types
    indrop_path = run_dir + '/' + INDROP_FILE
    if os.path.isfile(run_dir + '/' + INDROP_FILE):
            return 'indrop'
    return 'standard'

def check_sample_sheet(sample_sheet, run):
    # if sample sheet is not already there then copy the one from run_folder named
    # for flowcell
    if not os.path.exists(sample_sheet):
        flowcell = run.split('_')[-1][1:]
        path = config.SAMPLE_SHEET_DIR + flowcell + '.csv'
        if os.path.exists(path):
            util.copy(path, sample_sheet)

def write_cmd(cmd, output_dir, run):
    path = '%s/%s.opts' % (output_dir, run)
    with open(path, 'w') as fout:
        fout.write(cmd)

def bcl2fastq_process_runs():
    args, switches_to_names = initArgs()
    test = ('TEST' in args and args.TEST)
    run = os.path.basename(args.BCL_RUNFOLDER_DIR)
    setup_logging(run, test)
    test = False
    logging.info("***** START Odybcl2fastq *****\n\n")
    check_sample_sheet(args.BCL_SAMPLE_SHEET, run)
    logging.info("Beginning to process run: %s\n args: %s\n" % (run, json.dumps(vars(args))))
    sample_sheet = ss.sheet_parse(args.BCL_SAMPLE_SHEET)
    instrument = ss.get_instrument(sample_sheet['Data'])
    run_type = get_run_type(args.BCL_RUNFOLDER_DIR)
    mask_lists, mask_samples = extract_basemasks(sample_sheet['Data'], args.RUNINFO_XML, instrument, args, run_type)
    # skip everything but billing if run folder flagged
    if os.path.exists(args.BCL_RUNFOLDER_DIR + '/' + 'billing_only.txt'):
        logging.info("This run is flagged for billing only %s" % run)
        update_lims_db(run, sample_sheet, instrument)
        return
    jobs_tot = len(mask_lists)
    if jobs_tot > 1:
        logging.info("This run contains different masks in the same lane and will require %i bcl2fastq jobs" % jobs_tot)
    job_cnt = 1
    sample_sheet_dir = args.BCL_SAMPLE_SHEET
    output_dir = args.BCL_OUTPUT_DIR
    # run bcl2fatq per indexing strategy on run
    for mask, mask_list in mask_lists.items():
        output_suffix = None
        # if more than one bcl2fastq cmd needed suffix output dir and sample sheet
        if jobs_tot > 1:
            output_suffix = mask.replace(',', '_')
            args.BCL_OUTPUT_DIR = output_dir + '/' + output_suffix
            args.BCL_SAMPLE_SHEET = write_new_sample_sheet(mask_samples[mask], sample_sheet_dir, output_suffix)
        cmd = bcl2fastq_build_cmd(args,
                switches_to_names, mask_list, instrument, run_type)
        logging.info("\nJob %i of %i:" % (job_cnt, jobs_tot))
        if test:
            logging.info("Test run, command not run: %s" % cmd)
        else:
            logging.info('Launching bcl2fastq...%s\n' % cmd)
            success, message = bcl2fastq_runner(cmd,args)
            logging.info('message = %s' % message)
            summary_data = {}
            if success:
                # write bcl2fastq cmd
                write_cmd(cmd, args.BCL_OUTPUT_DIR, run)
                # update lims db
                update_lims_db(run, sample_sheet, instrument)
                # run  qc, TODO: consider a seperate job for this
                error_files, fastqc_err, fastqc_out = fastqc_runner(args.BCL_OUTPUT_DIR)
                output_log = get_output_log(run)
                with open(output_log, 'a+') as f:
                    f.write('\n'.join(fastqc_out) + "\n\n")
                    f.write('\n'.join(fastqc_err) + "\n\n")
                # copy run files to final
                copy_source_to_output(args.BCL_RUNFOLDER_DIR,
                        args.BCL_OUTPUT_DIR, args.BCL_SAMPLE_SHEET,
                        instrument)
                run_folder = args.BCL_OUTPUT_DIR.split('/').pop()
                # copy output to final dest where users will access
                copy_output_to_final(args.BCL_OUTPUT_DIR, run_folder, output_suffix)
                # get data from run to put in the email
                summary_data = parse_stats.get_summary(args.BCL_OUTPUT_DIR, instrument, args.BCL_SAMPLE_SHEET)
                summary_data['run'] = run
                summary_data['run_folder'] = run_folder
                summary_data['cmd'] = cmd
                summary_data['version'] = 'bcl2fastq2 v2.19'
                fastq_diff = compare_fastq(args.BCL_OUTPUT_DIR, instrument, run)
                logging.info('Fastq diff: %s' % json.dumps(fastq_diff))
            fromaddr = config.EMAIL['from_email']
            if success:
                toemaillist = config.EMAIL['to_email']
            else:
                toemaillist = config.EMAIL['to_email']
            logging.info('Sending email summary to %s\n' % json.dumps(toemaillist))
            sent = buildmessage(message, 'Demultiplex Summary for ' + run, summary_data, fromaddr, toemaillist)
            logging.info('Email sent: %s\n' % str(sent))
        job_cnt += 1
    logging.info("***** END Odybcl2fastq *****\n\n")
    return (0 if success else 1)

def get_output_log(run):
    return config.LOG_DIR + run + '.log'

def setup_logging(run, test):
    # take level from env or INFO
    level = os.getenv('LOGGING_LEVEL', logging.INFO)
    logging.basicConfig(
            filename=get_output_log(run),
            level=level,
            format='%(asctime)s %(filename)s %(message)s'
    )
    if test:
        logging.getLogger().addHandler(logging.StreamHandler())

if __name__ == "__main__":
    try:
        sys.exit(bcl2fastq_process_runs())
    except Exception as e:
        logging.exception(e)
        raise
