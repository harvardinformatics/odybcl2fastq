import unittest
import os
import csv
import json
from argparse import Namespace
import odybcl2fastq.util as util
import odybcl2fastq.run as run
import odybcl2fastq.parsers.parse_sample_sheet as ss

class Odybcl2fastqTests(unittest.TestCase):

    def setUp(self):
        self.sample_data_dir = (os.path.abspath( os.path.dirname( __file__ ) ) +
        '/sample_data/')

    def tearDown(self):
        pass

    def test_sheet_parse(self):
        sample_sheet_path = 'tests/sample_data/SampleSheet.csv'
        sample_sheet = ss.sheet_parse(sample_sheet_path)
        parts = ['Header', 'Reads', 'Settings', 'Data']
        for part in parts:
            assert (part in sample_sheet and sample_sheet[part])

    def test_get_instrument(self):
        run_info = 'tests/sample_data/RunInfo.xml'
        sample_sheet_path = 'tests/sample_data/SampleSheet.json'
        sample_sheet = util.load_json(sample_sheet_path)
        instrument =  ss.get_instrument(sample_sheet['Data'])
        assert instrument == 'hiseq'

    def test_extract_basemasks(self):
        run_info = 'tests/sample_data/RunInfo.xml'
        instrument = 'hiseq'
        sample_sheet_path = 'tests/sample_data/SampleSheet.json'
        sample_sheet = util.load_json(sample_sheet_path)
        mask_lists, mask_samples =  run.extract_basemasks(sample_sheet['Data'], run_info, instrument)
        mask_lists_control = {'y26,i8,y134': ['1:y26,i8,y134', '2:y26,i8,y134']}
        mask_samples_path = 'tests/sample_data/mask_samples.json'
        mask_samples_control = util.load_json(mask_samples_path)
        assert (mask_lists == mask_lists_control)
        assert (mask_samples == mask_samples_control)

    def test_build_cmd(self):
        mask_list = ['1:y26,i8,y134', '2:y26,i8,y134']
        instrument = 'hiseq'
        args = Namespace(BCL_ADAPTER_STRINGENCY=0.90000000000000002, BCL_BARCODE_MISMATCHES=0,
            BCL_CREATE_INDEXREAD_FASTQ=False, BCL_FASTQ_COMPRESSION_LEVEL=4,
            BCL_FIND_ADAPTERS_SLIDING_WINDOW=False, BCL_IGNORE_MISSING_BCLS=True,
            BCL_IGNORE_MISSING_FILTER=True, BCL_IGNORE_MISSING_POSITIONS=True,
            BCL_MASK_SHORT_ADAPTER_READS=22, BCL_MINIMUM_TRIMMED_READ_LENGTH=0,
            BCL_NO_BGZF=False, BCL_NO_LANE_SPLITTING=False,
            BCL_OUTPUT_DIR='/n/ngsdata/odybcl2fastq_test/171101_D00365_1013_AHYYTWBCXY',
            BCL_PROC_THREADS=8,
            BCL_RUNFOLDER_DIR='/n/boslfs/INSTRUMENTS/illumina/171101_D00365_1013_AHYYTWBCXY',
            BCL_SAMPLE_SHEET='/n/boslfs/INSTRUMENTS/illumina/171101_D00365_1013_AHYYTWBCXY/SampleSheet_new.csv',
            BCL_TILES=False, BCL_WITH_FAILED_READS=False,
            BCL_WRITE_FASTQ_REVCOMP=False,
            RUNINFO_XML='/n/boslfs/INSTRUMENTS/illumina/171101_D00365_1013_AHYYTWBCXY/RunInfo.xml',
            TEST=True
        )
        switches_to_names = {('--with-failed-reads',): 'BCL_WITH_FAILED_READS',
                ('--adapter-stringency',): 'BCL_ADAPTER_STRINGENCY', ('-p',
                    '--processing-threads'): 'BCL_PROC_THREADS', ('-o',
                        '--output-dir'): 'BCL_OUTPUT_DIR',
                    ('--find-adapters-with-sliding-window',):
                    'BCL_FIND_ADAPTERS_SLIDING_WINDOW',
                    ('--barcode-mismatches',): 'BCL_BARCODE_MISMATCHES',
                    ('--ignore-missing-positions',):
                    'BCL_IGNORE_MISSING_POSITIONS', ('--no-bgzf-compression',):
                    'BCL_NO_BGZF', ('--sample-sheet',): 'BCL_SAMPLE_SHEET',
                    ('--mask-short-adapter-reads',):
                    'BCL_MASK_SHORT_ADAPTER_READS',
                    ('--minimum-trimmed-read-length',):
                    'BCL_MINIMUM_TRIMMED_READ_LENGTH',
                    ('--ignore-missing-bcls',): 'BCL_IGNORE_MISSING_BCLS',
                    ('-R', '--runfolder-dir'): 'BCL_RUNFOLDER_DIR',
                    ('--create-fastq-for-index-reads',):
                    'BCL_CREATE_INDEXREAD_FASTQ',
                    ('--write-fastq-reverse-complement',):
                    'BCL_WRITE_FASTQ_REVCOMP', ('--no-lane-splitting',):
                    'BCL_NO_LANE_SPLITTING', ('--tiles',): 'BCL_TILES',
                    ('--ignore-missing-filter',): 'BCL_IGNORE_MISSING_FILTER',
                    ('--fastq-compression-level',):
                    'BCL_FASTQ_COMPRESSION_LEVEL'
        }
        run_type = None
        cmd_path = 'tests/sample_data/cmd.json'
        cmd_control = util.load_json(cmd_path)
        cmd = run.bcl2fastq_build_cmd(args,
                switches_to_names, mask_list, instrument, run_type)
        assert cmd == cmd_control

if __name__ == '__main__':
    unittest.main()
