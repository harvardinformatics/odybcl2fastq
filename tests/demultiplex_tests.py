import unittest
import os
import json
import gzip
import odybcl2fastq.parseargs as pa

class DemultiplexTests(unittest.TestCase):

    def setUp(self):
        self.sample_data_dir = (os.path.abspath( os.path.dirname( __file__ ) ) +
        '/sample_data/')

    def tearDown(self):
        pass

    def _load_json(self, path):
        obj = {}
        with open(path, 'r') as data:
            obj = json.load(data)
        return obj

    def test_demultiplex(self):
        cmd_path = 'tests/sample_data/cmd.json'
        cmd = self._load_json(cmd_path)
        #code, demult_out, demult_err = pa.run_cmd(cmd)
        #assert code == 0
        fastq_file = 'MDT1_SI_GA_A11_1_S1_L001_R1_001.fastq.gz'
        fastq_control_path = 'tests/sample_data/' + fastq_file
        #fastq_path = '/n/ngsdata/odybcl2fastq_test/171101_D00365_1013_AHYYTWBCXY/bambahmukku/' + fastq_file
        fastq_path = '/n/regal/informatics/mportermahoney/odytest/171101_D00365_1013_AHYYTWBCXY/bambahmukku/' + fastq_file
        control = gzip.open(fastq_control_path, 'r')
        test = gzip.open(fastq_path, 'r')
        for i in range(7):
            control_ln = next(control)
            test_ln = next(test)
            assert(control_ln == test_ln)

if __name__ == '__main__':
    unittest.main()
