from collections import OrderedDict
import logging
import odybcl2fastq.util as util
from odybcl2fastq.parsers.samplesheet import SampleSheet
from odybcl2fastq.parsers import parse_stats
import re, os

class Analysis(object):

    @staticmethod
    def create(run, mask, output_suffix):
        if run.type == '10x':
            return Analysis10x(run, mask, output_suffix)
        else:
            return AnalysisBcl2fastq(run, mask, output_suffix)

    def __init__(self, run, mask, output_suffix):
        self.run = run
        self.mask = mask
        if output_suffix:
            self.sample_sheet_path = self.run.sample_sheet.write_new_sample_sheet(self.mask, self.output_suffix)
        else:
            self.sample_sheet_path = self.run.sample_sheet_path
        self.sample_sheet = SampleSheet(self.sample_sheet_path)
        self.output_dir = run.output_dir
        if output_suffix:
            self.output_dir += '-' + output_suffix
        self.name = os.path.basename(self.output_dir)
        self.cmd = ''

    def run(self, output_log):
        # run unix cmd, stream out and error, return last lines of out
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

class Analysis10x(Analysis):

    def __init__(self, run, mask, output_suffix):
        super(AanalysisBcl2fastq, self).__init__()
        self.type = '10x'
        self.masklist = self.run.masklists[self.mask]

    def build_cmd(self, test = False):
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
                'run': self.run.args.BCL_RUNFOLDER_DIR,
                'samplesheet': self.sample_sheet_path,
                'output_dir': self.output_dir
                }
        cmd += ' '.join(['--' + k + '=' + v for (k, v) in options.items()])
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        path = self.write_cmd(cellranger, run.output_dir, run.name)
        slurm_cmd = 'sbatch %s' % path
        return slurm_cmd

class AnalysisBcl2fastq(Analysis):

    def __init__(self, run, mask, output_suffix):
        super(AnalysisBcl2fastq, self).__init__(run, mask, output_suffix)
        self.type = 'Bcl2fastq'
        self.version = 'bcl2fastq2 v2.2'
        self.mask_list = self.run.mask_lists[self.mask]

    def build_cmd(self, test = False):
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
        self.cmd = cmdstring
        self.write_cmd()
        return cmdstring



    def write_cmd(self):
        path = '%s/%s.opts' % (self.output_dir, self.name)
        with open(path, 'w') as fout:
            fout.write(self.cmd)

    def get_email_data(self):
        summary_data = parse_stats.get_summary(self.output_dir, self.run.instrument, self.sample_sheet_path, self.name)
        summary_data['run_folder'] = self.name
        summary_data['cmd'] = self.cmd
        summary_data['version'] = self.version




