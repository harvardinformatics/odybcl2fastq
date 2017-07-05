from lxml import html as lh
from collections import OrderedDict
import numpy

def get_lane_summary(path):
    data = parse_table(path)
    cols_to_display = ['lane', 'clusters', '% >= q30']
    filtered = []
    for r in data:
        row = OrderedDict()
        for k, v in r.items():
            if k in cols_to_display:
                row[k] = v.strip()
        filtered.append(row)
    return filtered

def parse_table(path):
    tree = lh.parse(path)
    table = tree.xpath("//h2[.='Lane Summary']/following::table[1]")[0]
    rows = iter(table)
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

def get_sample_summary(path):
    data = parse_table(path)
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
