import unittest
import os
from argparse import Namespace
import odybcl2fastq.run as run
from odybcl2fastq.parsers.samplesheet import SampleSheet
import logging

TEST_RUN_NAME       = 'sample_data'
BCL_OUTPUT_DIR      = 'test/%s/output_dir' % TEST_RUN_NAME
BCL_RUNFOLDER_DIR   = 'test/%s' % TEST_RUN_NAME
BCL_SAMPLE_SHEET    = 'test/%s/SampleSheet.csv' % TEST_RUN_NAME
RUNINFO_XML         = 'test/%s/RunInfo.xml' % TEST_RUN_NAME


class Odybcl2fastqTest(unittest.TestCase):

    def setUp(self):
        self.sample_data_dir = (
            os.path.abspath(os.path.dirname(__file__)) + '/sample_data/'
        )
        self.args = Namespace(
            BCL_ADAPTER_STRINGENCY=0.90000000000000002, BCL_BARCODE_MISMATCHES=0,
            BCL_CREATE_INDEXREAD_FASTQ=False, BCL_FASTQ_COMPRESSION_LEVEL=4,
            BCL_FIND_ADAPTERS_SLIDING_WINDOW=False, BCL_IGNORE_MISSING_BCLS=True,
            BCL_IGNORE_MISSING_FILTER=True, BCL_IGNORE_MISSING_POSITIONS=True,
            BCL_MASK_SHORT_ADAPTER_READS=22, BCL_MINIMUM_TRIMMED_READ_LENGTH=0,
            BCL_NO_BGZF=False, BCL_NO_LANE_SPLITTING=True,
            BCL_OUTPUT_DIR=BCL_OUTPUT_DIR,
            BCL_PROC_THREADS=8,
            BCL_RUNFOLDER_DIR=BCL_RUNFOLDER_DIR,
            BCL_SAMPLE_SHEET=BCL_SAMPLE_SHEET,
            BCL_TILES=False, BCL_WITH_FAILED_READS=False,
            BCL_WRITE_FASTQ_REVCOMP=False,
            RUNINFO_XML=RUNINFO_XML,
            TEST=True
        )
        self.switches_to_names = {('--with-failed-reads',): 'BCL_WITH_FAILED_READS',
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

    def tearDown(self):
        pass

    def testSheetParse(self):
        '''
        odybcl2fastq_tests: Test sample sheet parsing
        '''
        sample_sheet_path = 'test/sample_data/SampleSheet.csv'
        sample_sheet = SampleSheet(sample_sheet_path)
        parts = ['Header', 'Reads', 'Settings', 'Data']
        for part in parts:
            self.assertTrue(part in sample_sheet.sections and sample_sheet.sections[part])

    def testGetInstrument(self):
        '''
        odbcl2fastq_tests: Make sure the instrument can be retrieved from the sample sheet
        '''
        sample_sheet_path = 'test/sample_data/SampleSheet.csv'
        sample_sheet = SampleSheet(sample_sheet_path)
        instrument = sample_sheet.get_instrument()
        self.assertTrue(instrument == 'hiseq')

    def testExtractBasemasks(self):
        '''
        odbcl2fastq_tests: Make sure the base masks can be retrieved from the sample sheet
        '''
        run_info = 'test/sample_data/RunInfo.xml'
        instrument = 'hiseq'
        # json does not give ordered results
        sample_sheet_path = 'test/sample_data/SampleSheet.csv'
        sample_sheet = SampleSheet(sample_sheet_path)
        run_type = sample_sheet.get_run_type()
        mask_lists, mask_samples = run.extract_basemasks(sample_sheet.sections['Data'], run_info, instrument, self.args, run_type)
        mask_lists_control = {'y26,i8,y134': ['1:y26,i8,y134', '2:y26,i8,y134']}
        self.assertTrue(mask_lists == mask_lists_control)

    # def test_build_cmd(self):
    #     mask_list = ['1:y26,i8,y134', '2:y26,i8,y134']
    #     instrument = 'hiseq'
    #     run_type = None
    #     sample_sheet_path = 'tests/sample_data/SampleSheet.csv'
    #     sample_sheet = SampleSheet(sample_sheet_path)
    #     cmd_path = 'tests/sample_data/cmd.json'
    #     cmd_control = util.load_json(cmd_path)
    #     cmd = run.bcl2fastq_build_cmd(self.args,
    #             self.switches_to_names, mask_list, instrument, run_type, sample_sheet.sections)
    #     assert cmd == cmd_control

    def testProcessRuns(self):
        '''
        odybcl2fastq_tests: Use a dry run to make sure that logging output is correct
        '''
        logger = logging.getLogger('odybcl2fastq')
        run.bcl2fastq_process_runs(
            self.args,
            self.switches_to_names
        )
        self.assertTrue(os.path.exists(logger.handlers[0].baseFilename))
        odybcl2fastqlog = open(logger.handlers[0].baseFilename, 'r').read()
        self.assertTrue('Processing runs from run folder' in odybcl2fastqlog)

        runlogger = logging.getLogger('run_logger')
        self.assertTrue(os.path.exists(runlogger.handlers[0].baseFilename))
        runlog = open(runlogger.handlers[0].baseFilename, 'r').read()
        self.assertTrue('Beginning to process run' in runlog)
        self.assertTrue('END Odybcl2fastq' in runlog)


if __name__ == '__main__':
    unittest.main()
