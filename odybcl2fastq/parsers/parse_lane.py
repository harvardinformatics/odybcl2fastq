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
    for lane_info in data['ConversionResults']:
        lane = lane_info['LaneNumber']
        tot_yield = lane_info['Yield']
        clusters = lane_info['TotalClustersPF']
        reads = 0
        sam_num = 0
        lane_data = {'samples': [], 'clusters': clusters}
        # loop through samples and collect sam_data
        for row in lane_info['DemuxResults']:
            sam = row['SampleId']
            if sam not in summary_data:
                sam_num += 1
                sam_data = {
                        'sample': sam,
                        'indexes': [],
                        'yield': [],
                        'yield_q30': []
                }
                for i in row['IndexMetrics']:
                    sam_data['indexes'].append(i['IndexSequence'])
            else:
                sam_data = summary_data[sam]
            read = float(row['NumberReads'])
            reads += read
            sam_data['reads'] = read
            for r in row['ReadMetrics']:
                sam_data['yield'].append(float(r['Yield']))
                sam_data['yield_q30'].append(float(r['YieldQ30']))
            lane_data['samples'].append(sam_data)
        lane_data['reads'] = '{:,.0f}'.format(reads)
        lane_data['sam_num'] = sam_num
        summary_data[lane] = lane_data
    return summary_data

def format_lane_table(lanes):
    for lane_num, lane_info in lanes.items():
        for i, sam in enumerate(lane_info['samples']):
            row = OrderedDict()
            row['sample'] = sam['sample']
            row['index'] = ', '.join(sam['indexes'])
            row['reads'] = '{:,.0f}'.format(sam['reads'])
            row['% >= Q30'] = '{:.2f}'.format(numpy.sum(sam['yield_q30'])/numpy.sum(sam['yield']) * 100)
            lanes[lane_num]['samples'][i] = row
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
        for i, sam in enumerate(lane_info['samples']):
            sample = sam['sample']
            if sample not in agg:
                agg[sample] = {
                        'sample': sam['sample'],
                        'indexes': sam['indexes'],
                        'reads': 0,
                        'yield': [],
                        'yield_q30': []
                }
            agg[sample]['reads'] += sam['reads']
            agg_reads += sam['reads']
            agg[sample]['yield'].extend(sam['yield'])
            agg[sample]['yield_q30'].extend(sam['yield_q30'])
            lane_row['yield'].extend(sam['yield'])
            lane_row['yield_q30'].extend(sam['yield_q30'])
        lane_row['% >= Q30'] = numpy.sum(lane_row.pop('yield_q30'))/numpy.sum(lane_row.pop('yield'))
        lane_sum.append(lane_row)
    lanes_new = {}
    lanes_new[1] = lanes[1]
    lanes_new[1]['samples'] = agg.values()
    lanes_new[1]['reads'] = agg_reads
    ret = format_lane_table(lanes_new)
    ret[1]['lane_summary'] = lane_sum
    return ret
