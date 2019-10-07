#!python

#!/usr/bin/env python3
'''
Take from snakemake documentation
https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status

'''
import sys
import subprocess
if len(sys.argv) < 5:
    print("failed")
jobid = sys.argv[4]

cmd = "sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid
output = str(subprocess.check_output(cmd, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
failed_status=["FAILED", "OUT_OF_MEMORY", "TIMEOUT", "CANCELLED", "BOOT_FAIL", "DEADLINE", "NODE_FAIL", "PREEMPTED"]
if "COMPLETED" in output:
    print("success")
elif any(r in output for r in failed_status):
    print("failed")
else:
    print("running")
