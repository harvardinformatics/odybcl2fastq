from odybcl2fastq import util
import os

CONFIG_FILE = os.environ.get('ODY_CONFIG_FILE', os.path.join('config.json'))


class Config(object):

    def __init__(self):
        self.data = util.load_json(CONFIG_FILE)
        # add full paths for source, output and final for reusable scripts
        self.data['SOURCE_CLUSTER_PATH'] = '%s/' % os.environ.get('ODY_SOURCE')
        self.data['OUTPUT_CLUSTER_PATH'] = '%s/analysis/' % os.environ.get('ODY_SEQ_ROOT')
        self.data['PUBLISHED_CLUSTER_PATH'] = '%s/published/' % os.environ.get('ODY_SEQ_ROOT')
        self.data['REF_PATH'] = '%s/' % os.environ.get('ODY_REF')
        self.data['TEST'] = os.environ.get('ODY_TEST', 'FALSE') == 'TRUE'

    def __getattr__(self, attr):
        if attr in self.data:
            return self.data[attr]
        else:
            return None

    def __getitem__(self, key):
        return self.data[key]

    def keys(self):
        return self.data.keys()

    def __contains__(self, key):
        return key in self.data
