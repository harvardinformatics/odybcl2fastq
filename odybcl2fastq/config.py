from odybcl2fastq import util
import os

CONFIG_FILE = os.environ.get('ODYBCL2FASTQ_CONFIG_FILE', os.path.join('/etc', 'odybcl2fastq', 'config.json'))


class Config(object):

    def __init__(self):
        self.data = util.load_json(CONFIG_FILE)

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
