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
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from odybcl2fastq.parsers.makebasemask import extract_basemask

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
            'default'   : False,
            'action'    : 'store_true',
            'help'      : 'writes fastq files for index read(s)',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_BCLS',
            'switches'  : ['--ignore-missing-bcls'],
            'required'  : False,
            'help'      : 'missing or corrupt bcl files are ignored',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_FILTER',
            'switches'  : ['--ignore-missing-filter'],
            'required'  : False,
            'help'      : 'missing or corrupt filter files are ignored',
            'default'   : False,
            'action'    : 'store_true',
        },
        {
            'name'      : 'BCL_IGNORE_MISSING_POSITIONS',
            'switches'  : ['--ignore-missing-positions'],
            'required'  : False,
            'help'      : 'missing or corrupt positions files are ignored',
            'default'   : False,
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
        }
            
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
    
    for switches in switch_list:
        switch=[switch for switch in switches if '--' in switch][0]
        argvalue=str(argdict[switches_to_names[switches]])
        fout.write('%s\t%s\n' % (switch,argvalue))
        if argvalue !='False':
            cmdstrings.append(' '.join([switch,argvalue]))
   
    fout.close()
   
    cmdstring=' '.join(cmdstrings) 
   
    return cmdstring    
    



def main():
    bcl_namespace,attributedict,switches_to_names = initArgs()

    ### need to run check on basemask to see if can go with default pickup or need custom
    # newcmd will have to become a list resulting from iteration over make_bcl2fastq_cmd if
    # inconsistencies between runinfo.xml and sample_sheet.csv

    newcmd=make_bcl2fastq_cmd(attributedict,switches_to_names)
    print newcmd

if __name__ == "__main__":
    sys.exit(main())
