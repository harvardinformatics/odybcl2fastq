'''
odybcl2fastq Run

Created on  2018-12-20

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2018 The Presidents and Fellows of Harvard College.
All rights reserved.
@license: GPL v2.0
'''

import os
import hashlib
from glob import glob
import odybcl2fastq.util as util
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.parsers.makebasemask import extract_basemasks
from odybcl2fastq import config

class Run(object):
    """ Run object contains all command line args entered with details on what
    to run
    """
    def __init__(self, args, switches):
        self.args = args
        self.switches = switches
        self.test = ('TEST' in args and args.TEST)
        self.no_demultiplex = ('NO_DEMULTIPLEX' in args and args.NO_DEMULTIPLEX)
        self.no_post_process = ('NO_POST_PROCESS' in args and args.NO_POST_PROCESS)
        self.no_file_copy = ('NO_FILE_COPY' in args and args.NO_FILE_COPY)
        self.runfolder_dir = args.BCL_RUNFOLDER_DIR # full path
        self.name = os.path.basename(args.BCL_RUNFOLDER_DIR)
        self.sample_sheet_path = args.BCL_SAMPLE_SHEET
        self.output_dir = args.BCL_OUTPUT_DIR
        self.mask_lists = None
        self.mask_samples = None
        self.instrument = None
        self.type = None
        self.sample_sheet = None

    def set_sample_sheet(self):
        """ check for sample sheet existence, validate and correct, set
        instrument and type from sample sheet
        """
        # ensure sample sheet exists
        self.check_sample_sheet()
        self.sample_sheet = SampleSheet(self.sample_sheet_path)
        self.sample_sheet.validate()
        self.instrument = self.sample_sheet.get_instrument()
        # if output-dir was added to the sample sheet we need to set the
        # BCL_OUTPUT_DIR
        custom_output_dir = self.sample_sheet.get_output_dir()
        if custom_output_dir:
            self.output_dir = custom_output_dir
        self.type = self.sample_sheet.get_run_type()

    def check_sample_sheet(self):
        """ if doesn't exist, try to fetch by flowcell name
        """
        # if sample sheet is not already there then copy the one from run_folder named
        # for flowcell
        if not os.path.exists(self.sample_sheet_path):
            flowcell = self.name.split('_')[-1][1:]
            path = config.SAMPLE_SHEET_DIR + flowcell + '.csv'
            if os.path.exists(path):
                util.copy(path, self.sample_sheet_path)

    def fastq_checksum(self):
        """ get the checksum
        """
        checksum_report = self.output_dir + '/md5sum.txt'
        with open(checksum_report, 'a') as checksum_fh:
            checksum_fh.truncate(0)
            sample_proj_path = '%s/*/*.fastq.gz' % self.output_dir
            file_search = '%s/*.fastq.gz' % self.output_dir
            file_lst = glob(sample_proj_path)
            file_lst.extend(glob(file_search))
            for file_path in file_lst:
                with open(file_path, 'rb') as file_handle:
                    # put relative path in the checksum file so it can be run from
                    # the root even if the user copies the data somewhere else
                    checksum = hashlib.md5()
                    while True:
                        data = file_handle.read(2**20)
                        if not data:
                            break
                        checksum.update(data)
                    checksum = (checksum.hexdigest() + '  '
                                + file_path.replace(self.output_dir + '/', '') + '\n')
                    checksum_fh.write(checksum)

    def extract_basemasks(self):
        """ extract basemasks
        """
        self.mask_lists, self.mask_samples = extract_basemasks(
            self.sample_sheet.sections['Data'],
            self.args.RUNINFO_XML,
            self.instrument,
            self.args,
            self.type
        )
