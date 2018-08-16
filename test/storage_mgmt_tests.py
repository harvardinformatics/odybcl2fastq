import unittest
import subprocess
from odybcl2fastq import config
import datetime
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


class StorageMgmtTest(unittest.TestCase):

    def setUp(self):
        resetDirs(['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR'])

    def tearDown(self):
        resetDirs(['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR'])

    def testSetDirs(self):
        '''
        storage_mgmt_tests: Make sure the directories are properly set in the --help message
        '''
        cmd = 'storage_mgmt --help'
        code, out, err = run_cmd(cmd)
        self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
        storagepathkeys = ['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR']
        for k in storagepathkeys:
            self.assertTrue('(%s)' % config[k] in out, 'Cannot find %s in stdout: %s' % (config[k], out))

    def testMaxDelete(self):
        '''
        storage_mgmt_tests: Make sure that you don't delete more than max
        '''
        today = datetime.date.today()
        daysgoneby = 10
        rundate_str = (today - datetime.timedelta(days=daysgoneby)).strftime('%y%m%d')
        for n in xrange(1000, 1010):
            testrundir = '_'.join([rundate_str, 'D00742', str(n), 'BCBH04ANXX'])
            testrunpath = os.path.join(config.SOURCE_DIR, testrundir)
            os.makedirs(testrunpath)
        cmd = 'storage_mgmt --seq_storage SOURCE_DIR  --expired_after %d --delete' % (daysgoneby - 1)
        code, out, err = run_cmd(cmd)
        self.assertTrue(code > 0)
        self.assertTrue('No runs will be deleted.  You can override by passing in a higher max_delete parameter' in err, 'Incorrect output: %s' % err)

    def testLogging(self):
        '''
        storage_mgmt_tests: Make sure log file can be specified by env, either full path or relative to config LOG_DIR
        '''
        today = datetime.date.today()
        daysgoneby = 10
        rundate_str = (today - datetime.timedelta(days=daysgoneby)).strftime('%y%m%d')
        testrundir = '_'.join([rundate_str, 'D00742', '0197', 'BCBH04ANXX'])
        testrunpath = os.path.join(config.SOURCE_DIR, testrundir)
        os.makedirs(testrunpath)

        logfile = '/tmp/log'
        cmd = 'STORAGE_MGMT_LOG_FILE=%s storage_mgmt --seq_storage SOURCE_DIR  --expired_after %d' % (logfile, (daysgoneby - 1))
        code, out, err = run_cmd(cmd)
        self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
        logsdata = open(logfile, 'r').read()
        self.assertTrue('Found 1 runs for deletion' in logsdata, 'Incorrect output: %s' % logsdata)

        logfilename = 'test'
        logfile = os.path.join(config.LOG_DIR, logfilename)
        cmd = 'STORAGE_MGMT_LOG_FILE=%s storage_mgmt --seq_storage SOURCE_DIR  --expired_after %d' % (logfilename, (daysgoneby - 1))
        code, out, err = run_cmd(cmd)
        self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
        logsdata = open(logfile, 'r').read()
        self.assertTrue('Found 1 runs for deletion' in logsdata, 'Incorrect output: %s' % logsdata)

    def testNoDelete(self):
        '''
        storage_mgmt_tests: Test storage_mgmt dry run.  Should list directories to be removed.
        '''
        today = datetime.date.today()
        daysgoneby = 10
        rundate_str = (today - datetime.timedelta(days=daysgoneby)).strftime('%y%m%d')
        testrundir = '_'.join([rundate_str, 'D00742', '0197', 'BCBH04ANXX'])

        for d in ['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR']:
            testrunpath = os.path.join(config[d], testrundir)
            os.makedirs(testrunpath)

            cmd = 'storage_mgmt --seq_storage %s  --expired_after %d' % (d, (daysgoneby + 1))
            code, out, err = run_cmd(cmd)
            self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
            self.assertTrue('Found 0 runs for deletion' in err, 'Incorrect output: %s' % err)

            cmd = 'storage_mgmt --seq_storage %s  --expired_after %d' % (d, (daysgoneby - 1))
            code, out, err = run_cmd(cmd)
            self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
            self.assertTrue('Found 1 runs for deletion' in err, 'Incorrect output: %s' % err)

    def testDelete(self):
        '''
        storage_mgmt_tests: with --delete option on.  Should actually delete runs.
        '''
        today = datetime.date.today()
        daysgoneby = 10
        rundate_str = (today - datetime.timedelta(days=daysgoneby)).strftime('%y%m%d')
        testrundir = '_'.join([rundate_str, 'D00742', '0197', 'BCBH04ANXX'])
        for d in ['SOURCE_DIR', 'OUTPUT_DIR', 'FINAL_DIR']:
            testrunpath = os.path.join(config[d], testrundir)
            os.makedirs(testrunpath)

            cmd = 'storage_mgmt --seq_storage %s  --expired_after %d --delete' % (d, (daysgoneby + 1))
            code, out, err = run_cmd(cmd)
            self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
            self.assertTrue('Deleting 0 runs' in err, 'Incorrect output: %s' % err)

            cmd = 'storage_mgmt --seq_storage %s  --expired_after %d --delete' % (d, (daysgoneby - 1))
            code, out, err = run_cmd(cmd)
            self.assertTrue(code == 0, 'Error running command %s, %s' % (cmd, err))
            self.assertTrue('Deleting 1 runs' in err, 'Incorrect output: %s' % err)
            self.assertFalse(os.path.exists(testrunpath))
