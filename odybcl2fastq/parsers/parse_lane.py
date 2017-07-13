import os, re, glob
from lxml import html as lh
from collections import OrderedDict
import numpy

def get_summary(run_dir, short_id, instrument):
    if not os.path.exists(run_dir):
        raise UserException('Run directory does not exist: %s' % run_dir)
    summary_data = {'lanes':{}, 'reads':[]}
    for lane_dir in glob.glob(run_dir + "/Lane*/"):
        lane_data = {}
        lane_name =re.search('\/(Lane.*)\/', lane_dir).group(1)
        if instrument == 'nextseq':
            lane_file = lane_dir + 'html/' + short_id + '/all/all/all/lane.html'
            lane_data['reads'] = get_reads_nextseq(lane_file)
            lane_data['lane_summary'] = get_lane_summary_nextseq(lane_file)
            sample_file = lane_dir + 'html/' + short_id + '/all/all/all/laneBarcode.html'
            lane_data['sample_summary'] = get_sample_summary_nextseq(sample_file)
            lane_data['sample_num'] = get_sample_num(lane_data['sample_summary'])
        elif instrument == 'hiseq':
            lane_data['sample_summary'] = get_sample_summary_hiseq(path)
        else:
            raise Exception('instrument unknonw: ' + instrument)
        summary_data['lanes'][lane_name] = lane_data
    return summary_data

def get_sample_num(summary):
    sample_num = 0
    for row in summary:
        if row['index'].lower() not in ['unknown', 'undetermined']:
            sample_num += 1
    return sample_num

def get_reads_nextseq(path):
    data = parse_table_nextseq(path, 'Flowcell Summary')
    return data[0]['clusters(pf)']

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
    rows = iter(table)
    if 'Flowcell' not in tbl:
        next(rows) # skip first row which is high level header
    headers = [col.text.lower() for col in next(rows)]
    for i, h in enumerate(headers):
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
    for sam, data in sam_data.items():
        sam_data[sam]['reads'] = '{:,}'.format(int(numpy.sum(sam_data[sam]['reads'])))
        sam_data[sam]['% >= q30'] = '{:.2f}'.format(numpy.mean(sam_data[sam]['% >= q30']))
    return sam_data.values()
