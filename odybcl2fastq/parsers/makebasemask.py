from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException
from collections import OrderedDict,defaultdict
from copy import copy
from sets import Set

def make_universal_mask(rundata_by_read):
    universal_mask=OrderedDict()
    mask_type={'Y':'i','N':'y'}
    for read in rundata_by_read.keys():
        mask_value = rundata_by_read[read]['NumCycles']
        universal_mask[read] = '%s%s' % (mask_type[rundata_by_read[read]['IsIndexedRead']],mask_value)
    
    return universal_mask


def make_nextseq_mask(universal_mask,sample_dict):
    index1_length = int(universal_mask['read2'][1:])
    if len(sample_dict['index']) != index1_length:
        print('WARNING: sample 1st index of different length than determined from RunInfo.xml for %s' % sample_dict['Sample_Name'])        
        universal_mask['read2']='i%s' % len(sample_dict['index'])
    if 'i' in universal_mask['read3']:
        index2_length = int(universal_mask['read3'][1:])
        if len(sample_dict['index2']) != index2_length:
            universal_mask['read3']='i%s' % len(sample_dict['index2'])

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


def collapse_hiseq_masks(mask_lanes_dict,queues): # unique_masks is a python set
    """
    for hiseq runs, creates a list of lists, where each list represents a set
    of masks that can be run concurrently in a single demultiplexing run that
    obeys the following rule: each lane can have a separate mask but no lane
    can have multiple masks. this function build sets of masks for the minimum
    number of demultiplex runs necessary to handle all samples properly.
    """
 
    queue = {}
    masks = mask_lanes_dict.keys()
    for mask in masks:
        lanes = mask_lanes_dict[mask]
        for i,lane in enumerate(lanes):
            if lane not in queue.keys():
                queue[lane] = mask
                mask_lanes_dict[mask].pop(i)
        if len(mask_lanes_dict[mask]) == 0:
            mask_lanes_dict.pop(mask)                 
                    
    queues.append(queue) 
    if len(mask_lanes_dict) != 0:
        collapse_hiseq_masks(mask_lanes_dict,queues)
    
    queue_lists = []
    for queue in queues:
        queue_list = [key+':'+queue[key] for key in queue.keys()]
        queue_lists.append(queue_list)
    return queue_lists
    
     
def extract_basemasks(runinfo,sample_sheet):
    """
    creates a list of lists that contain masks
    that are compatible being run together in one demultiplexing
    instance, regardless of instrument; nextseq must run each different
    mask separately, while hiseq can handle different masks for different
    lanes but not different masks for the same lane
    """
    rundata_by_read = get_readinfo_from_runinfo(runinfo)
    data_by_sample = sheet_parse(sample_sheet)['Data']
    universal_mask=make_universal_mask(rundata_by_read)
    sample_masks=OrderedDict() # keep for figuring out which basemask runs to send to which clients ?
    unique_masks = Set()

    if 'Lane' in data_by_sample[data_by_sample.keys()[0]]:
        instrument = 'hiseq'
        for sample in data_by_sample.keys():
            lane,mask = make_hiseq_mask(universal_mask,sample,data_by_sample[sample])
            sample_masks[sample] = ':'.join([lane,mask])
            unique_masks.add(sample_masks[sample]) 
           
        if len(Set([mask.split(':')[1] for mask in unique_masks])) == 1:
            queue_lists = [list(list(Set([mask.split(':')[1] for mask in unique_masks])))]
        else:
            mask_lanes_dict = defaultdict(list)
            for mask in unique_masks:
                mask_lanes_dict[mask.split(':')[1]].append(mask.split(':')[0])
            queues = []
            queue_lists = collapse_hiseq_masks(mask_lanes_dict,queues)
                
    else:
        instrument = 'nextseq'
        for sample in data_by_sample.keys():
            sample_mask = make_nextseq_mask(universal_mask,data_by_sample[sample])    
            sample_masks[sample] = sample_mask
            unique_masks.add(sample_masks[sample])
        
            queue_lists = []
            for mask in unique_masks:
                queue_lists.append([mask])
                
                
    return queue_lists,instrument                
