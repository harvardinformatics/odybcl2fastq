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
import sys, os, re, traceback
import constants as const
from jinja2 import Environment, FileSystemLoader
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq.emailbuilder.emailbuilder import buildmessage
from odybcl2fastq.parsers import parse_lane
#emailer import buildmessage
from subprocess import Popen,PIPE
from os.path import basename

def initArgs():
    '''
    Setup arguments with parameterdef for bcl2fastq args
    other than the run and outpur dirs, check envs, parse
    commandline, return args
    '''

    parameterdefs = [
        {
            'name'      : 'BCL_LOAD_THREADS',
            'switches'  : ['-r','--loading-threads'],
            'required'  : False,
            'help'      : 'number threads to load bcl data',
            'default'   : 8,
            'type'      : int,
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
            'name'      : 'BCL_WRITE_THREADS',
            'switches'  : ['-w','--writing-threads'],
            'required'  : False,
            'help'      : 'number threads to write fastq files',
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
    bclargs = dict((attr, getattr(args,attr)) for attr in dir(args) if attr.startswith('BCL'))

    return args,bclargs,switches_to_names


def make_bcl2fastq_cmd(argdict,switches_to_names,runname='test'):
    cmdstrings=['bcl2fastq']
    fout=open('%s.opts' % runname,'w')
    # keeps consistent order of writing
    switch_list=switches_to_names.keys()
    switch_list.sort()
    #print 'argdict is', argdict
    for switches in switch_list:
        switch=[switch for switch in switches if '--' in switch][0]
        argvalue=str(argdict[switches_to_names[switches]])
        fout.write('%s\t%s\n' % (switch,argvalue))
        # the bit below prevents boolean flags from having values in the cmd
        if argvalue != 'False':
            if argvalue == 'True':
                cmdstrings.append(switch)
            else:
                cmdstrings.append(' '.join([switch,argvalue]))


    fout.close()

    cmdstring=' '.join(cmdstrings)

    return cmdstring


def bcl2fastq_build_cmd_by_queue():
    bcl_namespace,attributedict,switches_to_names = initArgs()
    newcmd=make_bcl2fastq_cmd(attributedict,switches_to_names)
    queuemasks,instrument =  extract_basemasks(bcl_namespace.RUNINFO_XML,bcl_namespace.BCL_SAMPLE_SHEET)
    cmds_by_queue = []
    for queue in queuemasks:
        if len(queuemasks) == 1:
            if instrument == 'nextseq':
                if len(queue) !=1:
                    raise UserException('more than 1 mask per queue detected for nextseq')

            queuecmd ='%s %s' % (newcmd, ' '.join(['--use-bases-mask %s' % mask for mask in queue]))

        elif instrument == 'hiseq':
            lanes = ','.join([mask.split(':')[0] for mask in queue])
            queuecmd ='%s --lanes %s %s' % (newcmd,lanes, ' '.join(['--use-bases-mask %s' % mask for mask in queue]))

        elif instrument == 'nextseq':
            if len(queue) == 1:
                queuecmd ='%s %s' % (newcmd,'--use-bases-mask %s' % queue[0])
            else:
                raise UserException('more than 1 mask per queue detected for nextseq')
        cmds_by_queue.append(queuecmd)
    return bcl_namespace,cmds_by_queue

def bcl2fastq_runner(cmd,bcl_namespace):
    demult_run = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    demult_out,demult_err=demult_run.communicate()
    if demult_run.returncode!=0:
        message = 'run %s failed\n%s\n' % (basename(bcl_namespace.BCL_RUNFOLDER_DIR),demult_err)
        print(message)
    else:
        message = 'run %s completed successfully\n' % basename(bcl_namespace.BCL_RUNFOLDER_DIR)

    fromaddr = 'adamfreedman@fas.harvard.edu'
    toemaillist=['adamfreedman@fas.harvard.edu']
    subject = basename(bcl_namespace.BCL_RUNFOLDER_DIR)
    lane_file = '/Users/portermahoney/tmp/analysis_finished/Lane1.indexlength_8/html/HVGC2BGX2/all/all/all/lane.html'
    lane_data = parse_lane.get_lane_summary(lane_file)
    sample_file = '/Users/portermahoney/tmp/analysis_finished/Lane1.indexlength_8/html/HVGC2BGX2/all/all/all/laneBarcode.html'
    sample_data = parse_lane.get_sample_summary(sample_file)
    # create html message with jinja
    j2_env = Environment(loader=FileSystemLoader(const.TEMPLATE_DIR),
            trim_blocks = True)
    # context to put in email template
    context = {
            'run': '170516_D00365_0940_BHFVGJBCXY',
            'subject': subject,
            'analyses': 1,
            'sample_no': len(sample_data),
            'reads': ['363','058','633'],
            'lane_data': lane_data,
            'sample_data': sample_data,
            'letter': '',
            # TODO: store in config or something
            'run_folder': 'http://software.rc.fas.harvard.edu/ngsdata/'

    }
    html = j2_env.get_template('summary.html').render(context)
    buildmessage(html,fromaddr,toemaillist,subject,ccemaillist=[],bccemaillist=[],server='smtp.fas.harvard.edu')

def bcl2fastq_process_runs(test=False):
    bcl_namespace,cmds = bcl2fastq_build_cmd_by_queue()
    for cmd in cmds:
        if test == True:
            print cmd
        else:
            print 'Launching bcl2fastq...%s\n' % cmd
            bcl2fastq_runner(cmd,bcl_namespace)


if __name__ == "__main__":
    #sys.exit(bcl2fastq_build_cmd_by_queue())
    sys.exit(bcl2fastq_process_runs())
