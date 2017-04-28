import xml.etree.ElementTree as ET

def get_readinfo_from_runinfo(runinfo_xml_file):
    tree = ET.parse(runinfo_xml_file)
    root = tree.getroot()

    read_dicts=[]

    for child in root.iter('Reads'):
        for child in child:
            read_dicts.append(child.attrib)

    return read_dicts

