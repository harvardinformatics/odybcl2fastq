from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo
from odybcl2fastq import UserException
from collections import OrderedDict,defaultdict
from copy import copy
import json
import re
import logging

NUTRAL_BASE = 'n'

logger = logging.getLogger('run_logger')

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
        index = re.sub('[^A-Z]+', '', index)
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

def make_mask(universal_mask, sample_key, sample_dict, log):
    # sample mask is based on universal mask from run info
    sample_mask = copy(universal_mask)
    # limit indexes to length in recipe
    if 'Recipe' in sample_dict and sample_dict['Recipe']:
        recipe = sample_dict['Recipe'].split('_')
        if 'index' in sample_dict:
            sample_dict['index'] = sample_dict['index'][:int(recipe[0])]
        if 'index2' in sample_dict and sample_dict['index2']:
            sample_dict['index2'] = sample_dict['index2'][:int(recipe[1])]
    # update index lengths from sample sheet
    if 'index' in sample_dict and sample_dict['index']: # both single and dual
        if 'read2' in sample_mask:
            sample_mask['read2'] = update_mask_index(sample_dict['index'],
                    sample_mask['read2'], sample_key)
    if 'index2' in sample_dict and sample_dict['index2'] and 'read3' in sample_mask: # dual indexed
        sample_mask['read3'] = update_mask_index(sample_dict['index2'],
                sample_mask['read3'], sample_key)
    elif 'read3' in sample_mask and 'i' in sample_mask['read3']:
        # if the universal mask had a indexed read3 but this sample doesn't, add nutral
        cnt = int(sample_mask['read3'][1:])
        sample_mask['read3'] = NUTRAL_BASE * cnt
    if log:
        logger.info('sample %s mask is: %s' % (sample_key, sample_mask))
    return ','.join(sample_mask.values())

def lists_from_mask(mask, data_by_sample):
    mask_lists = {}
    mask_samples = {}
    mask_lists[mask] = [mask]
    mask_samples[mask] = []
    for sample, row in data_by_sample.items():
        mask_samples[mask].append(row)
    return mask_lists, mask_samples

def extract_basemasks(data_by_sample, runinfo, instrument, run_type, log = True):
    """
    creates a list of lists that contain masks
    that are compatible being run together in one demultiplexing
    instance, regardless of instrument; nextseq must run each different
    mask separately, while hiseq can handle different masks for different
    lanes but if masks are different for same lane then it will require multiple
    demultiplexing runs
    """
    rundata_by_read = get_readinfo_from_runinfo(runinfo)
    universal_mask=make_universal_mask(rundata_by_read)
    mask_list = []
    mask_samples = {}
    mask_lists = {}
    if run_type == 'indrop':
        mask = []
        for read in universal_mask:
            mask.append('y*')
        mask = ','.join(mask)
        mask_lists, mask_samples = lists_from_mask(mask, data_by_sample)
    elif instrument in ['hiseq', 'novaseq']:
        # get mask per sample
        lane_masks=OrderedDict()
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row, log)
            lane = row['Lane']
            if lane not in lane_masks:
                lane_masks[lane] = set()
            lane_masks[lane].add(mask)
            if mask not in mask_samples:
                mask_samples[mask] = []
            mask_samples[mask].append(row)
            if log:
                logger.info('adding mask %s for lane %s' % (mask, lane))
        for lane, masks in lane_masks.items():
            for mask in list(masks):
                if mask not in mask_lists:
                    mask_lists[mask] = OrderedDict()
                mask_lists[mask][lane] = mask
        for mask, lane_masks in mask_lists.items():
            mask_lists[mask] = [lane + ':' + mask for lane, mask in lane_masks.items()]
    else: # nextseq miseq
        for sample, row in data_by_sample.items():
            mask = make_mask(universal_mask, sample, row, log)
            if mask not in mask_lists:
                mask_lists[mask] = [mask]
            if mask not in mask_samples:
                mask_samples[mask] = []
            mask_samples[mask].append(row)
    if log:
        logger.info('instrument is %s' % instrument)
        logger.info('masks: %s' % json.dumps(mask_lists))
    return mask_lists, mask_samples
