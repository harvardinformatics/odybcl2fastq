import xml.etree.ElementTree as ET
from collections import OrderedDict

def get_readinfo_from_runinfo(runinfo_xml_file):
    print('read')
    tree = ET.parse(runinfo_xml_file)
    root = tree.getroot()
    read_dicts=[]
    for child in root.findall('.//Read'):
        read_dicts.append(child.attrib)
    readkey_to_readdata_map=OrderedDict()
    for read_dict in read_dicts:
        number=read_dict.pop('Number')
        readkey_to_readdata_map['read%s' % number] = read_dict
    return readkey_to_readdata_map

