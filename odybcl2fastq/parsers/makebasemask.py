from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException
from collections import OrderedDict,defaultdict
from copy import copy
from sets import Set
import json

def make_universal_mask(rundata_by_read):
    universal_mask=OrderedDict()
    mask_type={'Y':'i','N':'y'}
    for read, data in rundata_by_read.items():
        mask_value = data['NumCycles']
        indexed = data['IsIndexedRead']
        if indexed not in mask_type:
            raise UserException('IsIndexedRead was a value not in %s: %s' % (json.dumps(mask_type), sample_dict['Sample_ID']))
        universal_mask[read] = '%s%s' % (mask_type[data['IsIndexedRead']],mask_value)
    return universal_mask


def make_nextseq_mask(universal_mask,sample_dict,legacytrim = True):
    #index1_length = int(universal_mask['read2'][1:])
    #if len(sample_dict['index']) != index1_length:
    #    print('WARNING: sample 1st index of different length than determined from RunInfo.xml for %s' % sample_dict['Sample_Name'])
    #    universal_mask['read2']='i%s' % len(sample_dict['index'])
    #if 'i' in universal_mask['read3']:
    #    index2_length = int(universal_mask['read3'][1:])
    #    if len(sample_dict['index2']) != index2_length:
    #        universal_mask['read3']='i%s' % len(sample_dict['index2'])
    universal_mask = copy(universal_mask)
    index1_length = int(universal_mask['read2'][1:])
    if len(sample_dict['index']) != index1_length:
        raise UserException('index 1 length inconsistent with sample: %s' % sample_dict['Sample_ID'])
    if 'i' in universal_mask['read3']:
        index2_length = int(universal_mask['read3'][1:])
        if len(sample_dict['index2']) != index2_length:
            raise UserException('index 2 length inconsistent with sample: %s' % sample_dict['Sample_ID'])
    if legacytrim == True:
        print 'universal mask is', universal_mask
        for read in universal_mask.keys():
            if 'y' == universal_mask[read][0]:
                newmask = 'y%sN' % str(int(universal_mask[read][1:])-1)
                universal_mask[read] = newmask
    return ','.join(universal_mask.values())

def update_mask_index(index, mask, sample_key):
    if 'i' in mask:
        bases = len(index)
        if bases <= 0:
            raise UserException('sample %s index has zero bases' % sample_key)
        mask = 'i%s' % bases
    else:
        raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample_key)
    return mask

def make_mask(universal_mask, sample_key, sample_dict):
    # sample mask is based on universal mask from run info
    sample_mask = copy(universal_mask)
    # update index lengths from sample sheet
    if 'index' in sample_dict: # both single and dual
        sample_mask['read2'] = update_mask_index(sample_dict['index'],
                sample_mask['read2'], sample_key)
    else:
        raise UserException('no index in sample sheet for %s' % sample_key)
    if 'index2' in sample_dict and sample_dict['index2']: # dual indexed
        print('dual index')
        sample_mask['read3'] = update_mask_index(sample_dict['index2'],
                sample_mask['read3'], sample_key)
    elif 'read3' in sample_mask and 'i' in sample_mask['read3']:
        # if the universal mask had a read3 but this sample doesn't, delete
        sample_mask['read3'] = sample_mask['read4']
        del(sample_mask['read4'])
    if sample_mask != universal_mask:
        print('sample mask for %s differs from the universal mask %s vs %s' % (sample_key, json.dumps(sample_mask), json.dumps(universal_mask)))
    return ','.join(sample_mask.values())

def make_hiseq_mask(universal_mask,sample_key,sample_dict,legacytrim = True):
    sample_mask = copy(universal_mask)
    print 'sampledict is', sample_dict
    if '_' not in sample_dict['Recipe']: # single indexed
        print('single')
        if len(sample_mask.keys()) in [2,3]: # single or paired end
            sample_mask['read2'] = 'i%s' % sample_dict['Recipe']
        else:
            raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample_key)
    elif '_' in sample_dict['Recipe']: # dual indexed
        print('dual')
        if 'i' in sample_mask['read2'] and 'i' in sample_mask['read3']:
            sample_mask['read2'] = 'i%s' % sample_dict['Recipe'].split('_')[0]
            sample_mask['read3'] = 'i%s' % sample_dict['Recipe'].split('_')[1]
        else:
            raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample_key)
    print(sample_mask)
    # TODO: bases 8 no, 6 yes add base for phasing
    # TODO: get rid of
    if legacytrim == True:
        for read in sample_mask.keys():
            if 'y' == sample_mask[read][0]:
                newmask = '%sN' % str(int(sample_mask[read][1:])-1)
            elif 'i' == sample_mask[read][0]:
                newmask = '%sN' % sample_mask[read]
                sample_mask[read] = newmask
            else:
                raise UserException('Unknown sample mask type for HiSeq sample %s' % sample_key)

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
    data_by_sample = sheet_parse(sample_sheet)['Data'] #list of tuples
    # TODO: universal_mask keep for figuring out which basemask runs to send to which clients ?
    universal_mask=make_universal_mask(rundata_by_read)
    if 'Lane' in data_by_sample[data_by_sample.keys()[0]]:
        instrument = 'hiseq'
        # get mask per sample
        lane_masks=OrderedDict()
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row)
            lane = row['Lane']
            if lane in lane_masks:
                # throw an error if mask is not the same per lane
                if lane_masks[lane] != mask:
                    raise UserException('hiseq sample mask for %s differs from another sample in lane %s: %s vs %s' % (sample_key, lane, mask, lane_mask[lane]))
            else:
                lane_masks[lane] = mask
                print('adding mask %s for lane %s' % (mask, lane))
        # one mask per lane
        mask_list = [lane + ':' + mask for lane, mask in lane_masks.items()]
    else:
        instrument = 'nextseq'
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row)
            if not mask_list:
                mask_list = [mask]
            elif mask not in mask_list:
                raise UserException('nextseq sample mask for %s differs from another sample %s vs %s' % (sample_key, mask, mask_list[0]))
    print('instrument is %s' % instrument)
    return mask_list, instrument
