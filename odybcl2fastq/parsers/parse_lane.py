import os, re, glob
from lxml import html as lh
from collections import OrderedDict
import constants as const
import numpy

def get_summary(run_dir, short_id, instrument, sample_sheet_dir):
    if not os.path.exists(run_dir):
        raise UserException('Run directory does not exist: %s' % run_dir)
    summary_data = {'lanes':{}}
    for lane_dir in glob.glob(run_dir + "/Lane*/"):
        lane_data = {}
        lane_name =re.search('\/(Lane.*)\/', lane_dir).group(1)
        if instrument == 'nextseq':
            lane_file = lane_dir + 'Reports/html/' + short_id + '/all/all/all/lane.html'
            lane_data['lane_summary'] = get_lane_summary_nextseq(lane_file)
            run = run_dir.split('/')[-1]
            sam_part = 'Reports/html/' + short_id + '/all/all/all/laneBarcode.html'
            sample_file = lane_dir + sam_part
            stats = run + '/' + lane_name + '/' + sam_part
            lane_data['sample_summary'], lane_data['reads'] = get_sample_summary_nextseq(sample_file)
            summary_data['undetermined'] = 'file labeled Undetermined_SO.'
        elif instrument == 'hiseq':
            stats = 'Basecall_Stats/Demultiplex_Stats.htm'
            sample_file = lane_dir + stats
            lane_data['sample_summary'], lane_data['reads'] = get_sample_summary_hiseq(sample_file)
            summary_data['undetermined'] = 'file(s) labeled Undetermined.'
        else:
            raise Exception('instrument unknonw: ' + instrument)
        lane_data['sample_num'] = get_sample_num(lane_data['sample_summary'])
        summary_data['stats_file'] = stats
        summary_data['sample_sheet'] = get_sample_sheet(sample_sheet_dir)
        summary_data['lanes'][lane_name] = lane_data
        summary_data['fastq_url'] = const.FASTQ_URL
        summary_data['fastq_dir'] = const.FASTQ_DIR
    return summary_data

def get_sample_sheet(sample_sheet_dir):
    data = ''
    with open(sample_sheet_dir, 'r') as ss:
        data = ss.read()
    return data

def get_sample_num(summary):
    sample_num = 0
    for row in summary:
        if row['index'].lower() not in ['unknown', 'undetermined']:
            sample_num += 1
    return sample_num

def get_lane_summary_nextseq(path):
    data = parse_table_nextseq(path, 'Lane Summary')
    cols_to_display = ['lane', 'clusters', '% >= q30']
    filtered = []
    for r in data:
        row = OrderedDict()
        for k, v in r.items():
            if k in cols_to_display:
                row[k] = v.strip()
        filtered.append(row)
    return filtered

def parse_table_nextseq(path, tbl):
    tree = lh.parse(path)
    table = tree.xpath("//h2[.='" + tbl + "']/following::table[1]")[0]
    table = table.xpath(".//tr") # get only rows
    rows = iter(table)
    next(rows) # skip first row which is high level header
    return rows_to_dict(rows)

def parse_table_hiseq(path, tbl):
    tree = lh.parse(path)
    # the html is hiseq file is malformed and headers are in seperate table from
    # data rows
    table = tree.xpath("//h2[.='" + tbl + "']/following::table[1]")[0]
    table2 = tree.xpath("//h2[.='" + tbl + "']/following::table[2]")[0]
    table = table.xpath(".//tr") # get only rows
    table2 = table2.xpath(".//tr") # get only rows
    rows = list(iter(table))
    rows.extend(list(iter(table2)))
    rows = iter(rows)
    return rows_to_dict(rows)

def rows_to_dict(rows):
    headers = [col.text.lower() for col in next(rows)]
    for i, h in enumerate(headers):
        # TODO: consider moving these to nextseq specific
        if h == '#':
            headers[i] = 'lane'
        if h == 'clusters':
            headers[i] = 'clusters raw'
            break
    data = []
    for r  in rows:
        values = [col.text for col in r]
        data.append(OrderedDict(zip(headers, values)))
    return data

def get_sample_summary_nextseq(path):
    data = parse_table_nextseq(path, 'Lane Summary')
    sam_data = OrderedDict()
    reads = 0
    for row in data:
        sam = row['sample']
        if sam not in sam_data:
            sam_data[sam] = OrderedDict()
            sam_data[sam]['sample'] = sam
            sam_data[sam]['index'] = row['barcode sequence']
            sam_data[sam]['reads'] = []
            sam_data[sam]['% >= q30'] = []
        sam_data[sam]['reads'].append(float(row['clusters'].replace(',', '')))
        sam_data[sam]['% >= q30'].append(float(row['% >= q30'].replace(',','')))
    reads = 0
    for sam, data in sam_data.items():
        read = int(numpy.sum(sam_data[sam]['reads']))
        reads += read
        sam_data[sam]['reads'] = '{:,.0f}'.format(read)
        sam_data[sam]['% >= q30'] = '{:.2f}'.format(numpy.mean(sam_data[sam]['% >= q30']))
    return sam_data.values(), '{:,.0f}'.format(reads)

def get_sample_summary_hiseq(path):
    data = parse_table_hiseq(path, 'Barcode lane statistics')
    sam_data = OrderedDict()
    reads = 0
    for row in data:
        sam = row['sample id']
        sam_data[sam] = OrderedDict()
        sam_data[sam]['sample'] = sam
        sam_data[sam]['index'] = row['index']
        read = float(row['# reads'].replace(',', ''))
        reads += read
        sam_data[sam]['reads'] = '{:,.0f}'.format(read)
        sam_data[sam]['% >= q30'] = '{:.2f}'.format(float(row['% of >= q30 bases (pf)'].replace(',','')))
    return sam_data.values(), '{:,.0f}'.format(reads)
