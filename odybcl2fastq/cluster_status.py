#!python

#!/usr/bin/env python3
'''
Take from snakemake documentation
https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status

'''
import sys
import subprocess

jobid = sys.argv[4]

cmd = "sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid
output = str(subprocess.check_output(cmd, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
    print("success")
elif any(r in output for r in running_status):
    print("running")
else:
    print("failed")
