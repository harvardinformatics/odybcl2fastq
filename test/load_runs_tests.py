import unittest
import subprocess
from odybcl2fastq import config
import shutil
import os


def run_cmd(cmd):
    # run unix cmd, return out and error
    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc.communicate()
    return (proc.returncode, out, err)


def resetDirs(dirs):
    for d in dirs:
        try:
            shutil.rmtree(config[d])
            os.makedirs(config[d])
        except Exception:
            pass


class LoadRunsTest(unittest.TestCase):

    def setUp(self):
        resetDirs(['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR'])

    def tearDown(self):
        resetDirs(['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR'])

    def testLogging(self):
        '''
        load_runs_tests: Make sure that log file specification works, either absolute path or relative path
        '''
        logfile = '/tmp/log'
        cmd = 'LOAD_RUNS_LOG_FILE=%s load_runs --no-daemon' % logfile
        code, out, err = run_cmd(cmd)
        self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
        logsdata = open(logfile, 'r').read()
        self.assertTrue('Found 0 need_to_process runs' in logsdata, 'Incorrect output: %s' % logsdata)
        self.assertTrue('Found 0 run_is_incomplete runs' in logsdata, 'Incorrect output: %s' % logsdata)
