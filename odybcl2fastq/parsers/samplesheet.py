from collections import OrderedDict
import logging
import odybcl2fastq.util as util
import re, os

class SampleSheet(object):

    def __init__(self, path):
        self.path = path
        self.lanes = []
        self.sections = self.sheet_parse(path)

    def sheet_parse(self, samplesheet=None):
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
        data_fields = []
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
                                lane = data_dict['Lane']
                            else:
                                name = '%s:%s' % (data_dict['Sample_Project'],data_dict['Sample_ID'])
                                lane = 1
                            if lane not in self.lanes:
                                self.lanes.append(lane)

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

    def get_instrument(self):
        if 'Lane' in self.sections['Data'].itervalues().next():
            instrument = 'hiseq'
        else:
            instrument = 'nextseq'
        return instrument

    def validate(self):
        if self.sample_names_corrected() or self.columns_corrected():
            # copy orig sample sheet as backup and record
            util.copy(self.path, self.path.replace('.csv', '_orig.csv'))
            # write a corrected sheet
            corrected_sample_sheet = self.write_new_sample_sheet(self.sections['Data'].values(), 'corrected')
            # copy corrected to sample sheet path, leave corrected file as record
            util.copy(corrected_sample_sheet, self.path)

    def columns_corrected(self):
        # remove all but whitelisted chars from data cols
        corrected = False
        invalid_regex = '[^\w\-@\. ]'
        for sam, line in self.sections['Data'].items():
            for col in line:
                if re.search(invalid_regex, line[col]):
                    corrected = True
                    tmp = line[col]
                    line[col] = re.sub(invalid_regex, '', line[col])
                    logging.info('Sample_Sheet corrected, invalid chars removed: %s to %s' % (tmp, line[col]))
            self.sections['Data'][sam] = line
        return corrected

    def sample_names_corrected(self):
        corrected = False
        proj_by_sample = {}
        cols_to_validate = ['Sample_ID', 'Sample_Name', 'Sample_Project']
        for sam, line in self.sections['Data'].items():
            for col in cols_to_validate:
                # replace any whitespace with underscores in sample names
                if util.contains_whitespace(line[col]):
                    corrected = True
                    tmp = line[col]
                    line[col] = line[col].replace(' ', '_')
                    logging.info('Sample_Sheet corrected, whitespace removed: %s to %s' % (tmp, line[col]))
                # remove any non alphanumeric chars
                if not util.alphanumeric(line[col]):
                    corrected = True
                    tmp = line[col]
                    line[col] = re.sub(r'[^\w-]', '', line[col])
                    if not line[col]:
                        raise Exception('For sample_sheet, %s: %s was all non alphanumeric' % (col, line[col]))
                    logging.info('Sample_Sheet corrected, non alphanumeric removed: %s to %s' % (tmp, line[col]))
            # if sample project is used, each sample id must belong to only one
            if line['Sample_Project']:
                if not line['Sample_ID'] in proj_by_sample:
                    proj_by_sample[line['Sample_ID']] = set()
                proj_by_sample[line['Sample_ID']].add(line['Sample_Project'])
            self.sections['Data'][sam] = line
        # rename samples that belong to multiple projects
        for id, projs in proj_by_sample.items():
            if len(projs) > 1:
                sam_proj_corrected = self.rename_samples(id)
                corrected = corrected and sam_proj_corrected
        return corrected

    def rename_samples(self, id):
        corrected = False
        # rename all samples with id by prefixing with submission
        for sam, line in self.sections['Data'].items():
            if line['Sample_ID'] == id:
                sub = line['Description']
                new_name = '%s_%s' % (sub, line['Sample_Name'])
                logging.info('renaming sample name: %s to %s' % (line['Sample_ID'],
                    new_name))
                line['Sample_ID'] = '%s_%s' % (sub, line['Sample_ID'])
                line['Sample_Name'] = new_name
                corrected = True
            self.sections['Data'][sam] = line
        return corrected

    def write_new_sample_sheet(self, new_samples, output_suffix):
        new_sample_sheet = self.path.replace('.csv', ('_' + output_suffix + '.csv'))
        input = open(self.path, 'r')
        output = open(new_sample_sheet, 'wb')
        for line in input:
            if not line.startswith('[Data]'):
                output.write(line)
            else:
                output.write(line)
                break
        # print data headers
        output.write(next(input))
        # write new samples to sheet
        new_lines = [(','.join(row.values()) + "\r\n") for row in new_samples]
        output.writelines(new_lines)
        output.close()
        input.close()
        return new_sample_sheet

    def get_output_dir(self):
        output_dir = ''
        if 'output-dir' in self.sections['Header']:
            output_dir = self.sections['Header']['output-dir'].strip()
        return output_dir

    def get_run_type(self):
        # chemistry field will be used for some special run types
        if 'Chemistry' in self.sections['Header'] and self.sections['Header']['Chemistry']:
            return self.sections['Header']['Chemistry'].strip().lower()
        else:
            return 'standard'

    def get_assay(self):
        if 'Assay' in self.sections['Header'] and self.sections['Header']['Assay']:
            return self.sections['Header']['Assay'].strip().lower()
        else:
            return None
