"""
function needs to be written so that an appropriate base mask 
can be written for each run/lane
"""
from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException
from collections import OrderedDict
from copy import copy
from sets import Set

def make_universal_mask(rundata_by_read):
    universal_mask=OrderedDict()
    mask_type={'Y':'i','N':'y'}
    for read in rundata_by_read.keys():
        mask_value = rundata_by_read[read]['NumCycles']
        universal_mask[read] = '%s%s' % (mask_type[rundata_by_read[read]['IsIndexedRead']],mask_value)
    
    return universal_mask


def make_nextseq_mask(universal_mask):
    return ','.join(universal_mask.values())

            
def make_hiseq_mask(universal_mask,sample_key,sample_dict):
    sample_mask = copy(universal_mask)
    
    if '_' not in sample_dict['Recipe']: # single indexed
        if len(sample_mask.keys()) in [2,3]: # single or paired end
            sample_mask['read2'] = 'i%s' % sample_dict['Recipe']
        else:
            raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample_key)
    elif '_' in sample_dict['Recipe']: # dual indexed
        if 'i' in sample_mask['read2'] and 'i' in sample_mask['read3']:
            sample_mask['read2'] = 'i%s' % sample_dict['Recipe'].split('_')[0]
            sample_mask['read3'] = 'i%s' % sample_dict['Recipe'].split('_')[1]
        else:
            raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample_key)            
    lane,mask = sample_dict['Lane'],','.join(sample_mask.values())

    return lane,mask            

"""
def collapse_hiseq_masks(sample_masks):
    number_lanes=int(sample_masks.keys()[-1].split(':'))
    masks_by_lane=OrderedDict()
    for mask in sample_masks.values():
        if mask.split(':')[1] not in masks_by_lane:
            masks_by_lane[mask.split(':')[1]]=Set()
        masks_by_lane[mask.split(':')[1]].add(mask.split(':')[0])

    if len(masks_by_lane) == 1: # only 1 mask for all lanes
        return masks_by_lane.keys()

    else: 
"""        
        
        
                

def extract_basemasks(runinfo,sample_sheet):
    rundata_by_read = get_readinfo_from_runinfo(runinfo)
    data_by_sample = sheet_parse(sample_sheet)['Data']
    universal_mask=make_universal_mask(rundata_by_read)
    sample_masks=OrderedDict() # keep for figuring out which basemask runs to send to which clients
    unique_masks = []

    if 'Lane' in data_by_sample[data_by_sample.keys()[0]]:
        unique_masks_by_lane = []
        for sample in data_by_sample.keys():
            lane,mask = make_hiseq_mask(universal_mask,sample,data_by_sample[sample])
            sample_masks[sample] = ':'.join([lane,mask])
            if sample_masks[sample] not in unique_masks_by_lane:
                unique_masks_by_lane.append(sample_masks[sample])            
                    
                                   
    else:
        for sample in data_by_sample.keys():
            sample_mask = make_nextseq_mask(universal_mask)    
            sample_masks[sample] = sample_mask
            if sample_masks[sample] not in unique_masks:
                unique_masks.append(sample_masks[sample])
        unique_masks_by_lane = unique_masks
                    
    return rundata_by_read,data_by_sample,universal_mask,sample_masks,unique_masks_by_lane        
  
   

