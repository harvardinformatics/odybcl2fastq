#!/usr/bin/env python3
'''
submit slurn jobs from snakemake, cp slurm cluster config into bash scripts

Created on  2019-03-25

@author: Meghan Correa <mportermahoney@g.harvard.edu>
@copyright: 2019 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

import os
import sys
import shutil
from collections import OrderedDict

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_props = read_job_properties(jobscript)

# uncomment out to write jobscript to script folder
#shutil.copyfile(jobscript, (job_props['input'][0] + '_jobscript'))

# the input file is a bash script to submit to slurm, read cmd in
with open(job_props['input'][0], 'r') as fh:
    cmd = fh.readlines()

prefix = []
cli_opts = []
# get slurm opts as a prefix for cmd
for prop, val in job_props['cluster'].items():
    sbatch_val = ''
    if val != '':
        sbatch_val = '=%s' % val
    prefix.append('#SBATCH --%s%s\n' % (prop, sbatch_val))
    cli_opts.append('--%s%s' % (prop, sbatch_val))

# write new prefixed cmd to the file
new_cmd = [cmd[0]] + prefix + cmd[1:]
with open(job_props['input'][0], 'w') as fh:
    for l in new_cmd:
        fh.write(l)


cli_opt_str = ' '.join(cli_opts)

os.system("sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript))
