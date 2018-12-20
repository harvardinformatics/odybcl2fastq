'''
odybcl2fastq Analysis

Created on  2018-12-20

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2018 The Presidents and Fellows of Harvard College.
All rights reserved.
@license: GPL v2.0
'''

import os
import sys
from subprocess import Popen
from subprocess import STDOUT
from collections import OrderedDict
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.parsers import parse_stats
from odybcl2fastq import config

class Analysis(object):
    """ analysis: suffix for 10x and for different basemasks, rewrites sample
    sheet with only samples for this analysis, building the command to be run
    """

    @staticmethod
    def create(run, mask, output_suffix):
        """ factory method to create analysis subclass
        """
        if run.type == '10x':
            # add 10x to the output suffix
            return Analysis10x(run, mask, output_suffix)
        return AnalysisBcl2fastq(run, mask, output_suffix)

    def __init__(self, run, mask, output_suffix):
        self.run = run
        self.output_dir = run.output_dir
        self.name = os.path.basename(self.output_dir)
        # suffix added to sample sheet and to output dir
        if output_suffix:
            sample_sheet_path = self.run.sample_sheet.write_new_sample_sheet(
                mask,
                output_suffix
            )
            self.output_dir += '-' + output_suffix
        else:
            sample_sheet_path = self.run.sample_sheet_path
        self.sample_sheet = SampleSheet(sample_sheet_path)
        self.cmd = ''
        self.template = None
        self.record_cmd = None

    def execute(self, output_log):
        """ execute command and stream output to log
        """
        # execute unix cmd, stream out and error, return last lines of out
        with open(output_log, 'a+') as writer:
            with open(output_log, 'r') as reader:
                proc = Popen(self.cmd, shell=True, stderr=STDOUT, stdout=writer)
                lines = []
                # stream output to log, std out
                while proc.poll() is None:
                    line = reader.readline()
                    # save last 40 lines for email
                    if line:
                        lines.append(line)
                        if len(lines) > 40:
                            lines.pop(0)
                    sys.stdout.write(line)
                # write any remaining output
                try:
                    line = reader.readline()
                    lines.append(line)
                    sys.stdout.write(line)
                except Exception:
                    pass
                code = proc.wait()
        return code, ''.join(lines)

    def get_template(self):
        """ get template
        """
        return self.template

class Analysis10x(Analysis):
    """ for 10x analysis
    """
    def __init__(self, run, mask, output_suffix):
        super(Analysis10x, self).__init__(run, mask, output_suffix)
        self.type = '10x'
        self.name += '_%s' % self.type
        self.template = 'summary10x.html'
        self.version = 'cellranger v2.1.0'

    def build_cmd(self, test=True):
        """ build 10x command for slurm and build slurm cmd to be run
            test defaulted to True for now
        """
        cmd = ("#!/bin/bash\n"
               "#SBATCH -p bos-info\n"
               "#SBATCH -N 1\n"
               "#SBATCH -n 8\n"
               "#SBATCH -t 08:00:00\n"
               "#SBATCH --mem 32000\n"
               "#SBATCH -J mkfastq\n"
               "#SBATCH -o cellranger_mkfastq_%A.out\n"
               "#SBATCH -e cellranger_mkfastq_%A.err\n")
        if test:
            cmd += "#SBATCH --test-only\n"

        cmd += ("source new-modules.sh\n"
                "module purge\n"
                "module load cellranger/2.1.0-fasrc01\n"
                "cellranger makfastq --localcores=8 --ignore-dual-index")
        options = {
            'run': self.run.runfolder_dir,
            'samplesheet': self.sample_sheet.path,
            'output_dir': self.output_dir
        }
        cmd += ' '.join(['--' + k + '=' + v for (k, v) in options.items()])
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.record_cmd = cmd
        path = self.write_cmd()
        slurm_cmd = 'sbatch %s' % path
        self.record_cmd += slurm_cmd
        return slurm_cmd

    def write_cmd(self):
        """ write cmd
        """
        path = '%s/%s_10x.sh' % (self.output_dir, self.name)
        with open(path, 'w') as fout:
            fout.write(self.record_cmd)
        return path

    def get_email_data(self):
        """ data for email template
        """
        summary_data = {}
        summary_data['run_folder'] = self.name
        summary_data['cmd'] = self.record_cmd
        summary_data['version'] = self.version
        summary_data['instrument'] = self.run.instrument
        summary_data['fastq_url'] = config.FASTQ_URL
        summary_data['sample_sheet'] = self.sample_sheet.path
        return summary_data


class AnalysisBcl2fastq(Analysis):
    """ for bcl2fastq analysis
    """
    def __init__(self, run, mask, output_suffix):
        super(AnalysisBcl2fastq, self).__init__(run, mask, output_suffix)
        self.type = 'Bcl2fastq'
        self.version = 'bcl2fastq2 v2.2'
        self.mask_list = self.run.mask_lists[mask]
        self.template = 'summary.html'

    def build_cmd(self, test=False):
        """ build bcl2fastq cmd
        """
        argdict = vars(self.run.args)
        cmdstrings = ['bcl2fastq']
        # keeps consistent order of writing
        switch_list = self.run.switches.keys()
        switch_list.sort()
        bcl_params = []
        cmd_dict = OrderedDict()
        for switches in switch_list:
            switch = [switch for switch in switches if '--' in switch][0]
            bcl_params.append(switch)
            argvalue = str(argdict[self.run.switches[switches]])
            # the bit below prevents boolean flags from having values in the cmd
            if argvalue != 'False':
                if argvalue == 'True':
                    cmd_dict[switch] = None
                else:
                    cmd_dict[switch] = argvalue
        if self.run.instrument == 'nextseq':
            cmd_dict['--no-lane-splitting'] = None
        # check for short reads, do not mask
        if self.run.type == 'indrop' or not self.sample_sheet.mask_short_reads():
            cmd_dict['--mask-short-adapter-reads'] = 0
        # grab any manually added params from sample sheet
        ss_params = self.sample_sheet.get_params_from_header(bcl_params)
        cmd_dict.update(ss_params)
        cmdstrings.extend([(k + ' ' + str(v)) if v is not None else k for k, v in cmd_dict.items()])
        # add the basemask generated from sample sheet if none was manually provided
        mask_switch = '--use-bases-mask'
        if mask_switch not in cmd_dict:
            # each mask should be prefaced by the switch
            mask_opt = mask_switch + ' ' + (' ' + mask_switch + ' ').join(self.mask_list)
            cmdstrings.append(mask_opt)
        cmdstring = ' '.join(cmdstrings)
        self.record_cmd = cmdstring
        self.write_cmd()
        return cmdstring

    def write_cmd(self):
        """ write cmd
        """
        path = '%s/%s.opts' % (self.output_dir, self.name)
        with open(path, 'w') as fout:
            fout.write(self.record_cmd)
        return path

    def get_email_data(self):
        """ data for email template
        """
        summary_data = parse_stats.get_summary(
            self.output_dir,
            self.run.instrument,
            self.sample_sheet.path,
            self.name
        )
        summary_data['run_folder'] = self.name
        summary_data['cmd'] = self.record_cmd
        summary_data['version'] = self.version
        return summary_data
