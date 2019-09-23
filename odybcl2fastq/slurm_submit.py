#!python

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

slurm_opts = OrderedDict([
        ('partition', '-p'),
        ('N', '-N'),
        ('cores', '-n'),
        ('time', '-t'),
        ('mem', '--mem'),
        ('name', '-J'),
        ('out', '-o'),
        ('error', '-e'),
        ('x', '-x')
])
jobscript = sys.argv[1]
job_props = read_job_properties(jobscript)
with open(job_props['input'][0] + 'test', 'w') as fh:
    with open(jobscript) as r:
        for l in r:
            if 'cd /app' in l or 'python3.6' in l:
                l = 'singularity exec -B /n/boslfs/INSTRUMENTS/illumina:/source \
                -B /n/boslfs/LABS/informatics/sequencing/PUBLISHED:/final \
                -B /n/boslfs/LABS/informatics/sequencing/ANALYSIS:/output \
                -B /n/informatics_external/seq/odybcl2fastq_log:/log \
                -B /n/boslfs/LABS/informatics/refs/10x/2019.05.19/cellranger:/ref \
                /n/boslfs/LABS/informatics/singularity_images/ody.sif ' + l
            fh.write(l)
shutil.copyfile((job_props['input'][0] + 'test'), jobscript)
# comment out the removal of the script to see the contents
os.remove((job_props['input'][0] + 'test'))

# the input file is a bash script to submit to slurm, read cmd in
with open(job_props['input'][0], 'r') as fh:
    cmd = fh.readlines()

prefix = []
cli_opts = []
# get slurm opts as a prefix for cmd
for prop, val in slurm_opts.items():
    if prop in job_props['cluster']:
        prefix.append('#SBATCH %s %s\n' % (val, job_props['cluster'][prop]))
        cli_opts.append('%s %s' % (val, job_props['cluster'][prop]))

# write new prefixed cmd to the file
new_cmd = [cmd[0]] + prefix + cmd[1:]
with open(job_props['input'][0], 'w') as fh:
    for l in new_cmd:
        fh.write(l)


cli_opt_str = ' '.join(cli_opts)

#os.system("hex exec sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript))
os.system("sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript))
