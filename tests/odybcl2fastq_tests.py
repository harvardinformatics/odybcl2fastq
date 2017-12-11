import unittest
import os
import csv
import json
import odybcl2fastq.parseargs as pa
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse

class Odybcl2fastqTests(unittest.TestCase):

    def setUp(self):
        self.sample_data_dir = (os.path.abspath( os.path.dirname( __file__ ) ) +
        '/sample_data/')

    def tearDown(self):
        pass

    def test_sheet_parse(self):
        sample_sheet_path = 'tests/sample_data/SampleSheet.csv'
        sample_sheet = sheet_parse(sample_sheet_path)
        parts = ['Header', 'Reads', 'Settings', 'Data']
        for part in parts:
            assert (part in sample_sheet and sample_sheet[part])

    def test_extract_basemasks(self):
        run_info = 'tests/sample_data/RunInfo.xml'
        sample_sheet_path = 'tests/sample_data/SampleSheet.json'
        sample_sheet = {}
        with open(sample_sheet_path, 'r') as data:
            sample_sheet = json.load(data)
        mask_lists, mask_samples, instrument =  pa.extract_basemasks(sample_sheet['Data'], run_info)
        assert instrument == 'hiseq'

    '''def test_build_cmd(self):
        cmd = pa.bcl2fastq_build_cmd(args,
                switches_to_names, mask_list, instrument, run_type)
        res = False
        with app.app_context():
            la = LipidAnalysis([self.sample_data_dir + 'neg_short.txt', self.sample_data_dir + 'pos_short.txt'])
            # check that rows have expected count
            res = (len(la.rows) == 64)
        assert res'''

if __name__ == '__main__':
    unittest.main()
