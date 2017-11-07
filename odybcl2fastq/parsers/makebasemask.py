from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq.parsers.parse_sample_sheet import sheet_parse
from odybcl2fastq import UserException
from collections import OrderedDict,defaultdict
from copy import copy
from sets import Set
import json
import logging

NUTRAL_BASE = 'n'

def make_universal_mask(rundata_by_read):
    universal_mask=OrderedDict()
    mask_type={'Y':'i','N':'y'}
    for read, data in rundata_by_read.items():
        mask_value = data['NumCycles']
        indexed = data['IsIndexedRead']
        if indexed not in mask_type:
            raise UserException('IsIndexedRead was a value not in %s: %s' % (json.dumps(mask_type), data['IsIndexedRead']))
        universal_mask[read] = '%s%s' % (mask_type[indexed],mask_value)
    return universal_mask

def update_mask_index(index, mask, sample):
    if 'i' in mask:
        bases = len(index)
        if bases <= 0:
            raise UserException('sample %s index has zero bases' % sample)
        prev_bases = int(mask.replace('i', ''))
        nutral_base = ''
        if prev_bases > bases:
            diff = prev_bases - bases
            nutral_base = NUTRAL_BASE * diff
        mask = 'i%s' % (str(bases) + nutral_base)
    else:
        raise UserException('# of runinfo reads inconsistent with sample sheet Recipe setting for %s' % sample)
    return mask

def make_mask(universal_mask, sample_key, sample_dict):
    # sample mask is based on universal mask from run info
    sample_mask = copy(universal_mask)
    # update index lengths from sample sheet
    if 'index' in sample_dict: # both single and dual
        if 'read2' in sample_mask:
            sample_mask['read2'] = update_mask_index(sample_dict['index'],
                    sample_mask['read2'], sample_key)
    else:
        raise UserException('no index in sample sheet for %s' % sample_key)
    if 'index2' in sample_dict and sample_dict['index2']: # dual indexed
        sample_mask['read3'] = update_mask_index(sample_dict['index2'],
                sample_mask['read3'], sample_key)
    elif 'read3' in sample_mask and 'i' in sample_mask['read3']:
        # if the universal mask had a indexed read3 but this sample doesn't, add nutral
        cnt = int(sample_mask['read3'][1:])
        sample_mask['read3'] = NUTRAL_BASE * cnt
    logging.info('sample %s mask is: %s' % (sample_key, sample_mask))
    if sample_mask != universal_mask:
        #TODO: consider removing this, too noisy
        logging.warning('sample mask for %s differs from the universal mask %s vs %s' % (sample_key, json.dumps(sample_mask), json.dumps(universal_mask)))
    return ','.join(sample_mask.values())

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
    universal_mask=make_universal_mask(rundata_by_read)
    mask_list = []
    if 'Lane' in data_by_sample.itervalues().next():
        instrument = 'hiseq'
        # get mask per sample
        lane_masks=OrderedDict()
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row)
            lane = row['Lane']
            if lane in lane_masks:
                # throw an error if mask is not the same per lane
                if lane_masks[lane] != mask:
                    raise UserException('hiseq sample mask for %s differs from another sample in lane %s: %s vs %s' % (sample, lane, mask, lane_masks[lane]))
            else:
                lane_masks[lane] = mask
                logging.info('adding mask %s for lane %s' % (mask, lane))
        # one mask per lane
        mask_list = [lane + ':' + mask for lane, mask in lane_masks.items()]
    else:
        instrument = 'nextseq'
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row)
            if not mask_list:
                mask_list = [mask]
            elif mask not in mask_list:
                raise UserException('nextseq sample mask for %s differs from another sample %s vs %s' % (sample, mask, mask_list[0]))
    logging.info('instrument is %s' % instrument)
    logging.info('mask is %s' % json.dumps(mask_list))
    return mask_list, instrument
