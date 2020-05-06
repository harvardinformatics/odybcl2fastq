#!/usr/bin/env python3

# -*- coding: utf-8 -*-

'''
Manage storage capacity by deleting old runs

Created on  2018-03-09
Updated on  2020-05-06

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2020 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import sys, os
import datetime
import json
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from subprocess import Popen, PIPE
from odybcl2fastq import config, initLogger


# Setup logging, to a file if set, otherwise, stderr
logger = initLogger('storage_mgmt', 'storage_mgmt')


# Set named storage paths
# If config file is identified, use those values
STORAGEPATHS = {
    'SOURCE_DIR' : '/source',
    'OUTPUT_DIR' : '/output',
    'PUBLISHED_DIR'  : '/published',
}
for k, v in STORAGEPATHS.items():
    STORAGEPATHS[k] = config[k]


def init_args():
    '''
    Setup arguments with parameterdef for bcl2fastq args
    other than the run and outpur dirs, check envs, parse
    commandline, return args
    '''

    parameterdefs = {
        'delete':
            {
                'required'  : False,
                'help'      : 'run in delete mode to actually delete runs',
                'action'    : 'store_true',
                'default'   : False,
            },
        'seq_storage':
            {
                'required'  : True,
                'choices'  : ['SOURCE_DIR', 'OUTPUT_DIR', 'PUBLISHED_DIR'],
                'help'      : 'Either the SOURCE_DIR ({SOURCE_DIR}), the OUTPUT_DIR ({OUTPUT_DIR}), or the PUBLISHED_DIR ({PUBLISHED_DIR})'.format(**STORAGEPATHS),
                'type'      : str,
            },
        'expired_after':
            {
                'required'  : True,
                'help'      : 'number of days after which runs can be deleted from the location given in seq_storage',
                'type'      : int,
            },
        'max_delete':
            {
                'required'  : False,
                'help'      : 'maximum number of runs to delete, if more runs are found then no runs will be deleted, this is a saftey measure',
                'type'      : int,
                'default'   : 5,
            }
    }

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

    # Use the parameterdefs for the ArgumentParser
    for switch, parameterdef in parameterdefs.items():
        if 'default' in parameterdef:
            parameterdef['help'] += '  [default: %s]' % parameterdef['default']
        parser.add_argument(*['--' + switch], **parameterdef)

    return parser.parse_args()


def run_cmd(cmd):
    # run unix cmd, return out and error
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    out, err = proc.communicate()
    return (proc.returncode, out, err)


def get_runs(seq_storage):
    # hard code paths as a safety measure so job cannot be used to delete
    # from unintended directories
    if seq_storage not in STORAGEPATHS:
        raise Exception('Unrecognized seq_storage option: %s' % seq_storage)

    path = STORAGEPATHS[seq_storage]
    logger.info('Cleaning up %s' % path)

    # get all the run folders
    cmd = 'find %s -maxdepth 1 -mindepth 1 -type d' % path
    code, out, err = run_cmd(cmd)
    if code != 0:
        raise Exception('Could not find directories at path %s: %s' % (cmd, err))
    return out.split()


def find_expired_runs(runs, oldest_str):
    to_delete = []
    for run in runs:
        run = run.decode('utf-8')
        run_parts = os.path.basename(run).split('_')
        run_date = run_parts[0]
        if len(run_parts) < 4 or len(run_date) != 6:
            logger.info('Badly formed run folder - ignoring %s' % run)
            continue
        if run_date < oldest_str:
            to_delete.append(run)
            logger.info('Adding %s to delete' % run)
    return to_delete


def manage_storage():
    args = init_args()
    runs = get_runs(args.seq_storage)

    # calculate expire date
    today = datetime.date.today()
    expire_str = (today - datetime.timedelta(days=args.expired_after)).strftime('%y%m%d')
    logger.info('Removing data older than %s' % expire_str)

    to_delete = find_expired_runs(runs, expire_str)
    delete_cnt = len(to_delete)

    # exit if we found more to delete than we planned as the max
    if delete_cnt > args.max_delete:
        logger.error('Runs to delete (%i) exceeds max_delete (%i).  No runs will be deleted.  You can override by passing in a higher max_delete parameter:\n%s' % (delete_cnt, args.max_delete, to_delete))
        sys.exit(1)

    # print out expired runs and delete if option is true
    to_delete_str = json.dumps(to_delete)
    cmds = []
    for dir in to_delete:
        cmds.append('rm -rf %s' % dir)
    logger.info('Delete cmds:\n%s' % json.dumps(cmds))
    if args.delete:
        logger.info('Deleting %i runs:\n%s' % (delete_cnt, to_delete_str))
        for cmd in cmds:
            code, out, err = run_cmd(cmd)
            if code != 0:
                raise Exception('Error deleting runs %s: %s' % (cmd, err))
        logger.info('Successfully deleted runs %s' % to_delete_str)
    else:
        logger.info('Found %i runs for deletion.  To delete pass in the parameter --delete:\n%s' % (delete_cnt, json.dumps(to_delete)))


if __name__ == "__main__":
    try:
        sys.exit(manage_storage())
    except Exception as e:
        logger.exception(e)
        sys.exit(1)
