import json
import logging
import requests
import os
from odybcl2fastq import config
from odybcl2fastq.parsers.parse_runinfoxml import get_readinfo_from_runinfo, get_runinfo
from odybcl2fastq.parsers.samplesheet import SampleSheet

class BauerDB(object):
    def __init__(self, sample_sheet_path):
        self.api = config.BAUER_DB['api']
        # use the seq app api for entering seq data
        self.seq_api = self.api + 'sequencing/api/'
        self.sample_sheet_path = sample_sheet_path
        self.token =  self.get_token()

    def insert_run(self):
        # insert run
        run = os.path.dirname(self.sample_sheet_path) + '/'
        runinfo_file = run + 'RunInfo.xml'
        run_data = get_runinfo(runinfo_file)
        run_id = self.post_data('runs', run_data)

        # insert read
        reads = get_readinfo_from_runinfo(runinfo_file)
        for i, read in enumerate(reads.values()):
            read_data = {
                    'run': run_id,
                    'number': i,
                    'indexed':  (1 if read['IsIndexedRead'] == 'Y' else 0),
                    'length': read['NumCycles']
            }
            read_id = self.post_data('reads', read_data)

        # insert lanes
        sample_sheet = SampleSheet(self.sample_sheet_path)
        lane_ids = {}
        for lane in sample_sheet.lanes:
            lane_data = {
                    'run': run_id,
                    'number': lane,
            }
            lane_ids[lane] = self.post_data('lanes', lane_data)

        # insert samples
        for sample_name, sample_row in sample_sheet.sections['Data'].items():
            sample_data = {
                    'name': sample_name,
                    'run': run_data['name'],
                    'description': sample_row['Description'],
                    'index1': sample_row['index'],
                    'index2': sample_row['index2']
            }
            # add a type if one was entered into sample sheet
            if 'Type' in sample_row and sample_row['Type']:
                sample_data['sample_type'] = self.get_sample_type(sample_row['Type'])
            if 'Lane' in sample_row and sample_row['Lane'].isdigit():
                lane = lane_ids[sample_row['Lane']]
            else:
                lane = lane_ids['1']
            sample_data['lane'] = lane
            sample_id = self.post_data('samples', sample_data)
        return True

    def get_token(self):
        return config.BAUER_DB['password']

    def post_data(self, endpoint, data):
        url = self.seq_api + endpoint + '/'
        print(url)
        headers = {'Authorization': 'Token %s' % self.token}
        r = requests.post(url = url, data = data, headers=headers)
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error('Api Error at %s: %s' % (url, e.response.text))
            raise e

        res = json.loads(r.text)
        item_id = res['id']
        logging.info('Loaded data for %s into id %d with data: %s' % (endpoint, item_id, json.dumps(data)))
        return item_id

    def get_data(self, endpoint):
        url = self.seq_api + endpoint + '/'
        headers = {'Authorization': 'Token %s' % self.token}
        r = requests.get(url = url, headers=headers)
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error('Api Error at %s: %s' % (url, e.response.text))
            raise e

        res = json.loads(r.text)
        logging.info('Retreived data for %s: %s' % (endpoint, json.dumps(res)))
        return res

    def get_sample_type(self, sample_type):
        sample_type = sample_type.lower()
        # get valid sample_types
        sample_types = self.get_data('sample_types')
        valid_sample_types = dict((t['name'].lower(), t['id']) for t in sample_types)
        # return valid sample_type or null
        if sample_type in valid_sample_types:
            return sample_type
        else:
            return None
