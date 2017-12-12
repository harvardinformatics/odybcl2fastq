from collections import OrderedDict
import logging

def sheet_parse(samplesheet=None):
    defaults_by_section = {
        'Header': {
        'IEMFileVersion': None,
        'InvestigatorName': None,
        'ExperimentName': None,
        'Date': None,
        'Workflow': None,
        'Application': None,
        'Assay': None,
        'Description': None,
        'Chemistry': None,},

        'Reads': {
        'read1_length': None,
        'read2_length': None},

        'Settings': {
        'Adapter': None,
        'AdapterRead2': None,
        'MaskAdapter': None,
        'MaskAdapterRead2': None,
        'FindAdaptersWithIndels': 1,
        'Read1EndWithCycle': None,
        'Read2EndWithCycle': None,
        'Read1StartFromCycle': None,
        'Read2StartFromCycle': None,
        'Read1UMILength': None,
        'Read2UMILength': None,
        'Read1UMIStartFromCycle': None,
        'Read2UMIStartFromCycle': None,
        'TrimUMI': 0,
        'ExcludeTiles': None,
        'CreateFastqForIndexReads': 0,
        'ReverseComplement': 0},

        'Data': OrderedDict(),
        }
    ssheet_open = open(samplesheet,'r')
    defaults_section = ''
    for line in ssheet_open:
        linelist = line.strip().split(',')
        if linelist[0] != '':
            if line[0] == '[':
                defaults_section = linelist[0][1:-1]
            else:
                if defaults_section in ['Settings','Header']:
                    defaults_by_section[defaults_section][linelist[0]] = linelist[1]

                elif defaults_section == 'Reads':
                    if defaults_by_section['Reads']['read1_length'] == None:
                        defaults_by_section['Reads']['read1_length'] = linelist[0]
                    elif linelist[0] != '':
                        defaults_by_section['Reads']['read2_length'] = linelist[0]

                else:
                    if 'Sample_ID' in linelist or 'SampleID' in linelist:
                        # TODO: lowercase all fields?
                        data_fields=[field.replace('SampleID','Sample_ID').replace('Index', 'index') for field in linelist if field != '']
                    else:
                        data_dict=OrderedDict(zip(data_fields,linelist[:len(data_fields)]))
                        if 'Lane' in data_dict.keys():
                            name = '%s:%s' % (data_dict['Lane'],data_dict['Sample_ID'])
                        else:
                            name = data_dict['Sample_ID']

                        defaults_by_section['Data'][name] = data_dict

    for section_key in ['Settings','Header','Reads']:
        for data_key in defaults_by_section[section_key].keys():
            if defaults_by_section[section_key][data_key] == None:
                defaults_by_section[section_key].pop(data_key)

    if len(defaults_by_section['Header']) == 0:
        #raise ValueError('No header information in sample sheet')
        logging.warning('odbcl2fastq WARNING: no header information in sample sheet')
    if len(defaults_by_section['Settings']) == 0:
        #raise ValueError('No settings information provided in sample sheet')
        logging.warning('odybcl2fastq WARNING: no settings information provided in sample sheet')

    if len(defaults_by_section['Reads']) == 0:
        #raise ValueError('No read information provided in sample sheet')
        logging.warning('odybcl2fastq WARNING: no read information provided in sample sheet')
    if len(defaults_by_section['Data']) == 0:
         raise ValueError('No data for samples present')
    return defaults_by_section

def get_instrument(sample_data):
    if 'Lane' in sample_data.itervalues().next():
        instrument = 'hiseq'
    else:
        instrument = 'nextseq'
    return instrument

