import os, re, glob
from lxml import html as lh
from odybcl2fastq import UserException
from collections import OrderedDict
import constants as const
import numpy
import json

# Stats.json
def get_summary(output_dir, short_id, instrument, sample_sheet_dir):
    stats_path = output_dir + '/Stats/Stats.json'
    if not os.path.exists(stats_path):
        raise UserException('Stats path does not exist: %s' % stats_path)
    with open(stats_path) as f:
        data = json.load(f)
    if not data:
        raise UserException('Stats file empty: %s' % stats_path)
    summary_data = {} # used in email
    lanes = get_stats(data)
    if instrument == 'hiseq':
        lanes = format_lane_table(lanes)
    elif instrument == 'nextseq':
        lanes = format_nextseq_tables(lanes)
    else:
        raise Exception('instrument unknonw: ' + instrument)
    summary_data['lanes'] = lanes
    summary_data['stats_file'] = stats_path
    summary_data['sample_sheet'] = get_sample_sheet(sample_sheet_dir)
    summary_data['fastq_url'] = const.FASTQ_URL
    summary_data['fastq_dir'] = const.FASTQ_DIR
    return summary_data

def get_sample_sheet(sample_sheet_dir):
    data = ''
    # putting sample sheet in the email is a convenience so don't fail if path
    # for sample sheet is wrong
    if os.path.exists(sample_sheet_dir):
        with open(sample_sheet_dir, 'r') as ss:
            data = ss.read()
    return data

def get_stats(data):
    summary_data = OrderedDict()
    # get data on from each lane and samples in the lane
    for lane_info in data['ConversionResults']:
        lane = lane_info['LaneNumber']
        sam_num = 0
        lane_stats = {
                'samples': OrderedDict(),
                'clusters': lane_info['TotalClustersPF'],
                'reads': lane_info['Yield'],
                'sam_num': len(lane_info['DemuxResults'])
        }
        # loop through samples and collect sam_data
        for row in lane_info['DemuxResults']:
            sam = row['SampleId']
            lane_stats['samples'][sam] = get_sam_stats(sam, lane_stats['samples'], row)
        # include a row for undertermined stats
        if lane_info['Undetermined']:
            sam = 'undetermined'
            lane_stats['samples'][sam] = get_sam_stats(sam, lane_stats['samples'], lane_info['Undetermined'])
        summary_data[lane] = lane_stats
    return summary_data

def get_sam_stats(sam, sam_summary, row):
    if sam not in sam_summary:
        sam_stats = {
                'sample': sam,
                'indexes': [],
                'yield': [],
                'yield_q30': []
        }
        if 'IndexMetrics' in row:
            for i in row['IndexMetrics']:
                sam_stats['indexes'].append(i['IndexSequence'])
        else:
            sam_stats['indexes'].append('undetermined')
    else:
        sam_stats = sam_summary[sam]
    sam_stats['reads'] = float(row['NumberReads'])
    for r in row['ReadMetrics']:
        sam_stats['yield'].append(float(r['Yield']))
        sam_stats['yield_q30'].append(float(r['YieldQ30']))
    return sam_stats

def format_lane_table(lanes):
    for lane_num, lane_info in lanes.items():
        for sam_name, sam_info in lane_info['samples'].items():
            row = OrderedDict()
            row['sample'] = sam_info['sample']
            row['index'] = ', '.join(sam_info['indexes'])
            row['reads'] = '{:,.0f}'.format(sam_info['reads'])
            row['% >= Q30'] = '{:.2f}'.format(numpy.sum(sam_info['yield_q30'])/numpy.sum(sam_info['yield']) * 100)
            lanes[lane_num]['samples'][sam_name] = row
    return lanes

def format_nextseq_tables(lanes):
    lane_sum = []
    agg = OrderedDict()
    agg_reads = 0
    for lane_num, lane_info in lanes.items():
        lane_row = OrderedDict()
        lane_row['lane'] = lane_num
        lane_row['clusters'] = lane_info['clusters']
        lane_row['yield'] = []
        lane_row['yield_q30'] = []
        for sam_name, sam_info in lane_info['samples'].items():
            if sam_name not in agg:
                agg[sam_name] = {
                        'sample': sam_info['sample'],
                        'indexes': sam_info['indexes'],
                        'reads': 0,
                        'yield': [],
                        'yield_q30': []
                }
            agg[sam_name]['reads'] += sam_info['reads']
            agg_reads += sam_info['reads']
            agg[sam_name]['yield'].extend(sam_info['yield'])
            agg[sam_name]['yield_q30'].extend(sam_info['yield_q30'])
            lane_row['yield'].extend(sam_info['yield'])
            lane_row['yield_q30'].extend(sam_info['yield_q30'])
        lane_row['% >= Q30'] = numpy.sum(lane_row.pop('yield_q30'))/numpy.sum(lane_row.pop('yield'))
        lane_sum.append(lane_row)
    lanes_new = {}
    lanes_new[1] = lanes[1]
    lanes_new[1]['samples'] = agg
    lanes_new[1]['reads'] = agg_reads
    ret = format_lane_table(lanes_new)
    ret[1]['lane_summary'] = lane_sum
    return ret
