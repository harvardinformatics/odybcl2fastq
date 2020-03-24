from collections import OrderedDict
import logging
import odybcl2fastq.util as util
import pandas as pd
import re, os

class SampleSheet(object):
    SAMPLE_SHEET_FILE = 'SampleSheet.csv'

    def __init__(self, path):
        self.path = path
        self.lanes = []
        self.sections = self.sheet_parse(path)
        self.samples = self.get_samples()

    def get_samples(self):
        data = pd.DataFrame(self.sections['Data'].values())
        return data

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
                            if 'Lane' in list(data_dict):
                                name = '%s:%s' % (data_dict['Lane'],data_dict['Sample_ID'])
                                lane = data_dict['Lane']
                            else:
                                name = '%s:%s' % (data_dict['Sample_Project'],data_dict['Sample_ID'])
                                lane = '1'
                            if lane not in self.lanes:
                                self.lanes.append(lane)

                            defaults_by_section['Data'][name] = data_dict

        for section_key in ['Settings','Header','Reads']:
            for data_key in list(defaults_by_section[section_key]):
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
        run = os.path.basename(os.path.dirname(self.path))
        instrument_name = run.split('_')[1]
        instrument = ''
        if instrument_name.startswith('D'):
            instrument = 'hiseq'
        elif instrument_name.startswith('N'):
            instrument = 'nextseq'
        elif instrument_name.startswith('A'):
            instrument = 'novaseq'
        elif instrument_name.startswith('M'):
            instrument = 'miseq'
        else:
            raise ValueError('Instrument %s does not match known types: D, N, A, M' % instrument_name)
        return instrument

    def validate(self):
        if self.validate_sample_names() | self.validate_index2():
            # copy orig sample sheet as backup and record
            util.copy(self.path, self.path.replace('.csv', '_orig.csv'))
            # write a corrected sheet
            corrected_sample_sheet = self.write_new_sample_sheet(self.sections['Data'].values(), 'corrected')
            # copy corrected to sample sheet path, leave corrected file as record
            util.copy(corrected_sample_sheet, self.path)

    def validate_index2(self):
        corrected = False
        fixed = self.samples
        if 'index2' in fixed.columns:
            if fixed[fixed.index2!=''].empty:
                # we have a totally empty index2 column so delete
                fixed = fixed.drop('index2', axis = 1)
                if 'I5_index_ID' in fixed.columns:
                    fixed = fixed.drop('I5_index_ID', axis = 1)
                corrected = True
                fixed_rows = fixed.to_dict('records')
                for i, k in enumerate(self.sections['Data'].keys()):
                    self.sections['Data'][k] = OrderedDict(fixed_rows[i])
        return corrected

    def validate_sample_names(self):
        corrected = False
        proj_by_sample = {}
        for sam, line in self.sections['Data'].items():
            cols_to_validate = ['Sample_ID', 'Sample_Name', 'Sample_Project']
            for col in cols_to_validate:
                if col in line:
                    # remove any whitespace
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
        new_sample_sheet = os.path.dirname(self.path) + '/' + self.SAMPLE_SHEET_FILE
        new_sample_sheet = new_sample_sheet.replace('.csv', ('_' + output_suffix + '.csv'))
        input = open(self.path, 'r')
        output = open(new_sample_sheet, 'w')
        for line in input:
            if not line.startswith('[Data]'):
                output.write(line)
            else:
                output.write(line)
                break
        # write new samples to sheet
        new_lines = []
        for i, row in enumerate(new_samples):
            # print data headers
            if i == 0:
                new_lines.append(','.join(row.keys()) + "\r\n")
            new_lines.append(','.join(row.values()) + "\r\n")
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

    def get_sample_types(self):
        types = {}
        for key, row in self.sections['Data'].items():
            name = key.split(':')[1]
            if 'Type' in row and row['Type']:
                types[name] = row['Type']
        return types

    def get_sample_projects(self):
        sp = {}
        for key, row in self.sections['Data'].items():
            name = key.split(':')[1]
            if 'Sample_Project' in row and row['Sample_Project']:
                if row['Sample_Project'] not in sp:
                    sp[row['Sample_Project']] = []
                sp[row['Sample_Project']].append(name)
        return sp

    def get_assay(self):
        if 'Assay' in self.sections['Header'] and self.sections['Header']['Assay']:
            return self.sections['Header']['Assay'].strip().lower()
        else:
            return None

    def get_submissions(self):
        instrument = self.get_instrument()
        subs = set()
        for key, row in self.sections['Data'].items():
            if row['Description']:
                subs.add(row['Description'])
        return list(subs)

    def has_poly_A_index(self):
        for k, row in self.sections['Data'].items():
            if 'AAAA' in row['index']:
                return True
            else:
                return False
