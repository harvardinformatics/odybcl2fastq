#!/usr/bin/env python

# -*- coding: utf-8 -*-

'''
manage storage capacity by deleting old runs

Created on  2018-03-09

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2018 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''
import sys, os
import datetime
import logging
import json
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from collections import OrderedDict
import odybcl2fastq.util as util
from odybcl2fastq import constants as const
from odybcl2fastq import config
from subprocess import Popen, PIPE, STDOUT

LOG_FILE = const.ROOT_DIR + 'storage_mgmt.log'

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
        'path':
            {
                'required'  : True,
                'help'      : 'location to check storage',
                'type'      : str,
            },
        'expired_after':
            {
                'required'  : True,
                'help'      : 'number of days after which runs can be deleted from the location given in path',
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
        parser.add_argument(*['--' + switch],**parameterdef)

    return parser.parse_args()

def run_cmd(cmd):
    # run unix cmd, return out and error
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    out, err = proc.communicate()
    return (proc.returncode, out, err)


def setup_logging():
    # take level from env or INFO
    level = os.getenv('LOGGING_LEVEL', logging.INFO)
    logging.basicConfig(
            filename=LOG_FILE,
            level=level,
            format='%(asctime)s %(filename)s %(message)s'
    )
    # add to standard out for testing purposes
    logging.getLogger().addHandler(logging.StreamHandler())

def get_runs(path):
    # get all the run folders
    cmd = 'find %s -maxdepth 1 -mindepth 1 -type d' % path
    code, out, err = run_cmd(cmd)
    if code != 0:
        raise Exception('Could not fine directories at path %s: %s' % (cmd, err))
    return out.split()

def find_expired_runs(runs, oldest_str):
    to_delete =  []
    for run in runs:
        run_parts = os.path.basename(run).split('_')
        run_date = run_parts[0]
        if len(run_parts) < 4 or len(run_date) != 6:
            logging.info('Baddly formed run folder - ignoring %s' % run)
            continue
        if run_date < oldest_str:
            to_delete.append(run)
            logging.info('Adding %s to delete' % run)
    return to_delete

def manage_storage():
    setup_logging()
    args = init_args()
    runs = get_runs(args.path)

    # calculate expire date
    today = datetime.date.today()
    expire_str = (today - datetime.timedelta(days=args.expired_after)).strftime('%y%m%d')

    to_delete = find_expired_runs(runs, expire_str)
    delete_cnt = len(to_delete)

    # exit if we found more to delete than we planned as the max
    if delete_cnt > args.max_delete:
        raise Exception('Runs to delete (%i) exceeds max_delete (%i).  No runs will be deleted.  You can override by passing in a higher max_delete parameter' % (delete_cnt, args.max_delete))

    # print out expired runs and delete if option is true
    if args.delete:
        to_delete_str = json.dumps(to_delete)
        logging.info('Deleting %i runs:\n%s' % (delete_cnt, to_delete_str))
        cmd = 'rm -rf {%s}' % ','.join(to_delete)
        code, out, err = run_cmd(cmd)
        if code != 0:
            raise Exception('Error deleting runs %s: %s' % (cmd, err))
        logging.info('Successfully deleted runs %s' % to_delete_str)
    else:
        logging.info('Found %i runs for deletion.  To delete pass in the parameter --delete:\n%s' % (delete_cnt, json.dumps(to_delete)))

if __name__ == "__main__":
    try:
        sys.exit(manage_storage())
    except Exception as e:
        logging.exception(e)
        raise
