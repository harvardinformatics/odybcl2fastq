import os, re
from odybcl2fastq import UserException
from collections import OrderedDict
from odybcl2fastq import config
import operator
import locale
import numpy
import json
import logging

MIN_UNDETER_CNT = 1000000

def get_summary(output_dir, instrument, sample_sheet_dir, run_folder):
    """
    parse summary from Stats.json
    """
    # set locale so numeric formating works
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    stats_url = config.FASTQ_URL + run_folder + '/Reports/html'
    stats_path = output_dir + '/Stats/Stats.json'
    if not os.path.exists(stats_path):
        raise UserException('Stats path does not exist: %s' % stats_path)
    with open(stats_path) as f:
        data = json.load(f)
    if not data:
        raise UserException('Stats file empty: %s' % stats_path)
    # parse stats for both hiseq and nextseq
    lanes = get_stats(data)
    lane_sum = []
    if instrument in ['nextseq', 'miseq']:
        # no real lanes in nextseq, aggregate the stats
        lane_sum = get_lane_sum(lanes)
        lanes = aggregate_nextseq_lanes(lanes)
        undetermined = format_undetermined_nextseq(data['UnknownBarcodes'])
    elif instrument == 'hiseq':
        undetermined = format_undetermined(data['UnknownBarcodes'])
    else:
        raise Exception('instrument unknonw: ' + instrument)
    # format lane summary tables
    lanes = format_lane_table(lanes)
    summary_data = {
            'lane_sum': lane_sum,
            'lanes': lanes,
            'instrument': instrument,
            'stats_file': stats_url,
            'sample_sheet': get_sample_sheet(sample_sheet_dir),
            'fastq_url': config.FASTQ_URL,
            'fastq_dir': config.PUBLISHED_CLUSTER_PATH,
            'undetermined': undetermined,
            'undetermined_file': 'Undetermined_SO',
            'sample_sheet_file': sample_sheet_dir
    }
    logging.info("summary_data for email: %s\n" % json.dumps(summary_data))
    return summary_data

def get_stats(data):
    stats = OrderedDict()
    # get data on from each lane and samples in the lane
    for lane_info in data['ConversionResults']:
        lane = lane_info['LaneNumber']
        lane_stats = {
                'samples': OrderedDict(),
                'clusters': lane_info['TotalClustersPF'],
                'yield': lane_info['Yield'],
                'sam_num': len(lane_info['DemuxResults'])
        }
        # loop through samples and collect sam_data
        for row in lane_info['DemuxResults']:
            sam = row['SampleId']
            lane_stats['samples'][sam] = get_sam_stats(sam, lane_stats['samples'], row)
        # include a row for undertermined stats
        if 'Undetermined' in lane_info:
            sam = 'undetermined'
            lane_stats['samples'][sam] = get_sam_stats(sam, lane_stats['samples'], lane_info['Undetermined'])
        stats[lane] = lane_stats
    return stats

def get_sam_stats(sam, sam_summary, row):
    if sam not in sam_summary:
        sam_stats = {
                'sample': sam,
                'index': [],
                'yield': [],
                'yieldq30': []
        }
        if 'IndexMetrics' in row:
            for i in row['IndexMetrics']:
                sam_stats['index'].append(i['IndexSequence'])
        else:
            sam_stats['index'].append('undetermined')
    else:
        sam_stats = sam_summary[sam]
    sam_stats['reads'] = float(row['NumberReads'])
    for r in row['ReadMetrics']:
        # ignore 0s
        if float(r['Yield']) > 0.0:
            sam_stats['yield'].append(float(r['Yield']))
            sam_stats['yieldq30'].append(float(r['YieldQ30']))
    return sam_stats

def get_lane_sum(lanes):
    lane_sum = []
    for lane_num, info in lanes.items():
        reads = 0
        lane_yield = []
        lane_yieldq30 = []
        for sam in info['samples'].values():
            reads += sam['reads']
            lane_yield.extend(sam['yield'])
            lane_yieldq30.extend(sam['yieldq30'])
        q30 = numpy.sum(lane_yieldq30) / numpy.sum(lane_yield) * 100
        lane_row = OrderedDict()
        lane_row['lane'] = lane_num
        lane_row['clusters'] = locale.format('%d', reads, True)
        lane_row['% Bases >= Q30'] = locale.format('%.2f',  q30, True)
        lane_sum.append(lane_row)
    return lane_sum

def format_lane_table(lanes):
    for lane_num, lane_info in lanes.items():
        lanes[lane_num]['clusters'] = locale.format('%d', lane_info['clusters'], True)
        for sam_name, sam_info in lane_info['samples'].items():
            row = OrderedDict()
            row['sample'] = sam_info['sample']
            row['index'] = ', '.join(sam_info['index'])
            row['clusters'] = locale.format('%d', sam_info['reads'], True)
            q30 = numpy.sum(sam_info['yieldq30']) / numpy.sum(sam_info['yield']) * 100
            row['% >= Q30'] = locale.format('%.2f', q30, True)
            lanes[lane_num]['samples'][sam_name] = row
    return lanes

def aggregate_nextseq_lanes(lanes):
    lanes_new = {}
    agg = OrderedDict()
    agg_reads = 0
    # aggregate sample stats since there are no real lanes in nextseq
    for lane_num, lane_info in lanes.items():
        for sam_name, sam_info in lane_info['samples'].items():
            if sam_name not in agg:
                agg[sam_name] = {
                        'sample': sam_info['sample'],
                        'index': sam_info['index'],
                        'reads': 0,
                        'yield': [],
                        'yieldq30': []
                }
            agg[sam_name]['reads'] += sam_info['reads']
            agg_reads += sam_info['reads']
            agg[sam_name]['yield'].extend(sam_info['yield'])
            agg[sam_name]['yieldq30'].extend(sam_info['yieldq30'])
    # since nextseq has no lanes put aggregated results in a single "lane"
    lanes_new[1] = {
            'sam_num': lanes[1]['sam_num'],
            'samples': agg,
            'clusters': agg_reads
    }
    return lanes_new

def format_undetermined_nextseq(undeter):
    all = {}
    top = OrderedDict()
    for lane in undeter:
        for index, cnt in lane['Barcodes'].items():
            if index not in all:
                all[index] = 0
            all[index] += cnt
    sorted_undeter = OrderedDict(sorted(all.items(), key=lambda kv: kv[1], reverse = True))
    for index, cnt in sorted_undeter.items():
        if cnt > MIN_UNDETER_CNT:
                top[index] = locale.format('%d', cnt, True)
        else:
            break
    return {1: top}

def format_undetermined(undeter):
    top = OrderedDict()
    for lane in undeter:
        top[lane['Lane']] = OrderedDict()
        sorted_undeter = sorted(lane['Barcodes'].items(), key=operator.itemgetter(1), reverse = True)
        for (index, cnt) in sorted_undeter:
            if cnt > MIN_UNDETER_CNT:
                top[lane['Lane']][index] = locale.format('%d', cnt, True)
    return top
