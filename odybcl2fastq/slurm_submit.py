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
        ('chdir', '--chdir')
])
jobscript = sys.argv[1]
job_props = read_job_properties(jobscript)
with open(job_props['input'][0] + 'test', 'w') as fh:
    with open(jobscript) as r:
        for l in r:
            if 'cd /app' in l or 'python3.6' in l:
                l = 'singularity exec -B /n/boslfs/INSTRUMENTS/illumina:/source \
                -B /n/boslfs/LABS/informatics/sequencing/PUBLISHED/odybcl2fastq_test:/final \
                -B /n/boslfs/LABS/informatics/sequencing/ANALYSIS/odytest:/output \
                -B ~/odybcl2fastq_log_test:/log \
                -B /n/boslfs/LABS/informatics/refs/10x/2019.05.19/cellranger:/ref \
                -B /etc/slurm:/etc/slurm -B /slurm:/slurm \
                -B /usr/bin/sbatch:/usr/bin/sbatch -B /usr/bin/sacct:/usr/bin/sacct \
                -B /n/informatics/repos/odybcl2fastq_10x/odybcl2fastq:/app \
                -B /etc/sssd/sssd.conf:/etc/sssd/sssd.conf \
                -B /usr/lib64/slurm:/usr/lib64/slurm \
                -B /usr/lib64/libmunge.so.2:/usr/lib64/libmunge.so.2 \
                -B /usr/lib64/libmunge.so.2.0.0:/usr/lib64/libmunge.so.2.0.0 \
                -B /usr/lib64/libslurmdb.so.33:/usr/lib64/libslurmdb.so.33 \
                -B /var/run/munge:/var/run/munge \
                -B /etc/nsswitch.conf:/etc/nsswitch.conf \
                -B /etc/sssd/:/etc/sssd/ \
                -B /usr/share/lua:/usr/share/lua \
                -B /var/lib/sss:/var/lib/sss ~/ody_dev.sif ' + l
            fh.write(l)
shutil.copyfile((job_props['input'][0] + 'test'), jobscript)

# the input file is a bash script to submit to slurm, read cmd in
with open(job_props['input'][0], 'r') as fh:
    cmd = fh.readlines()

prefix = []
cli_opts = []
# get slurm opts as a prefix for cmd
for prop, val in slurm_opts.items():
    prefix.append('#SBATCH %s %s\n' % (val, job_props['cluster'][prop]))
    cli_opts.append('%s %s' % (val, job_props['cluster'][prop]))

# write new prefixed cmd to the file
new_cmd = [cmd[0]] + prefix + cmd[1:]
with open(job_props['input'][0], 'w') as fh:
    for l in new_cmd:
        fh.write(l)


cli_opt_str = ' '.join(cli_opts)

#os.system("hex exec sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript))
str_cmd = "sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript)
with open(job_props['input'][0] + 'cmd', 'w') as fh:
    fh.write(str_cmd)
os.system("sbatch {cli_opts} {script}".format(cli_opts=cli_opt_str, script=jobscript))
