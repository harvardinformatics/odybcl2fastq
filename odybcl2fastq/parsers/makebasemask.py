"""
function needs to be written so that an appropriate base mask 
can be written for each run/lane
"""

from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException

def check_runinfo_samplesheet_consistency(readdict,sheetdict):
    index_lengths={'index1':'','index2':''}
    for read in readdict.keys():
        if read == 'read2' and readict[read]['IsIndexedRead'] == 'Y':
            index_lengths['index1'] = int(readdict[read]['NumCycles'])
        elif read == 'read3' and readict[read]['IsIndexedRead'] == 'Y':
            index_lengths['index2'] = int(readdict[read]['NumCycles'])     

    for sample in sheetdict['Data'].keys():
        if sheetdict['Data'][sample]['index'] == '' and index_lengths['index1']!= '':
            raise UserException('Index sequence not present in %s' % sample)
        elif len(sheetdict['Data'][sample]['index']) != index_lengths['index1']:
            return False
        elif sheetdict['Data'][sample]['index2'] == '' and index_lengths['index2']!= '':
            raise UserException('Index2 sequence not present in %s' % sample)
        elif len(sheetdict['Data'][sample]['index2']) != index_lengths['index2']:
            return False
        else return True

def extract_basemask(runinfo,sample_sheet):
    rundata_by_read=get_readinfo_from_runinfo(runinfo)
    data_by_sample=sheet_parse(sample_sheet)['Data']
    
    masktest=check_runinfo_samplesheet_consistency(rundata_by_read,data_by_sample_sheet)
    if masktest == True:
        return None
    elif masktest == False:
        return False
        
        

 
