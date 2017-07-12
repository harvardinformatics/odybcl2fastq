'''
Created on Jul 1, 2015

@author: afreedman
'''
import unittest
import os
import sys
import pickle
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse

data_dir='test_samplesheets/'

class Test(unittest.TestCase):
		
    def setUp(self):
        self.segood_hiseq_pickledict = '1_hiseq_singleindex_pe_samplesheet.dict.pickled' 
        self.segood_hiseq_ssheet = '1_hiseq_singleindex_pe_samplesheet.csv'
        self.pegood_singleindex_nextseq_ssheet='2_SampleSheet_NextSeq_PE_singleindex_170518_NS500422_0501_AHVFNVBGX2.csv'
        self.pegood_singleindex_nextseq_pickledict='2_SampleSheet_NextSeq_PE_singleindex_170518_NS500422_0501_AHVFNVBGX2.dict.pickled'
    def tearDown(self):
        pass

    def testValidateSampleSheetParse(self):
        picklefile='%s%s' % (data_dir,self.segood_hiseq_pickledict)
        benchmark = pickle.load(open(picklefile,'rb'))
        testparse = sheet_parse('%s%s' % (data_dir,self.segood_hiseq_ssheet)) 
        self.assertEqual(benchmark,testparse)

    def testValidateNextSeqPESingleIndexGood(self):
        picklefile='%s%s' % (data_dir,self.pegood_singleindex_nextseq_pickledict)
        benchmark2 = pickle.load(open(picklefile,'rb'))
        testparse2 = sheet_parse('%s%s' % (data_dir,self.pegood_singleindex_nextseq_ssheet))
        self.assertEqual(benchmark2,testparse2)

if __name__ == "__main__":
    unittest.main()
