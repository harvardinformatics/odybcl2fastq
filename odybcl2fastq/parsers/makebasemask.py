"""
function needs to be written so that an appropriate base mask 
can be written for each run/lane
"""

from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException

def check_runinfo_samplesheet_consistency(readdict,sampledict):
    index_lengths={'index1':'','index2':''}
    for read in readdict.keys():
        if read == 'read2' and readdict[read]['IsIndexedRead'] == 'Y':
            index_lengths['index1'] = int(readdict[read]['NumCycles'])
        elif read == 'read3' and readdict[read]['IsIndexedRead'] == 'Y':
            index_lengths['index2'] = int(readdict[read]['NumCycles'])     
    
    print 
    for sample in sampledict.keys():
        if sampledict[sample]['index'] == '' and index_lengths['index1']!= '':
            raise UserException('Index sequence not present in %s' % sample)
        elif len(sampledict[sample]['index']) != index_lengths['index1']:
            return False
        elif sampledict[sample]['index2'] == '' and index_lengths['index2']!= '':
            raise UserException('Index2 sequence not present in %s' % sample)
        elif len(sampledict[sample]['index2']) != index_lengths['index2']:
            return False
        else: 
            return True

def extract_basemask(runinfo,sample_sheet):
    rundata_by_read=get_readinfo_from_runinfo(runinfo)
    data_by_sample=sheet_parse(sample_sheet)['Data']
    
    masktest=check_runinfo_samplesheet_consistency(rundata_by_read,data_by_sample)
    if masktest == True:
        return True
    elif masktest == False:
        return False
        
        

 
