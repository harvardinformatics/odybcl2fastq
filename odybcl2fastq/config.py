from odybcl2fastq import util
import os

CONFIG_FILE = os.environ.get('ODY_CONFIG_FILE', os.path.join('config.json'))
DEFAULT_INFORMATICS_ROOT = '/n/boslfs02/LABS/informatics'
DEFAULT_SEQ_ROOT = '%s/sequencing' % DEFAULT_INFORMATICS_ROOT


class Config(object):

    def __init__(self):
        self.data = util.load_json(CONFIG_FILE)
        # add full paths for source, output and final for reusable scripts
        self.data['SEQ_ROOT'] = self.check_dir('%s/' % os.environ.get('ODY_SEQ_ROOT', '%s/' % DEFAULT_SEQ_ROOT).rstrip('/'))
        self.data['SOURCE_CLUSTER_PATH'] = self.check_dir('%s/' % os.environ.get('ODY_SOURCE', '%s/source' % DEFAULT_SEQ_ROOT).rstrip('/'))
        self.data['OUTPUT_CLUSTER_PATH'] = self.check_dir('%s/analysis/' % os.environ.get('ODY_SEQ_ROOT', DEFAULT_SEQ_ROOT).rstrip('/'))
        self.data['PUBLISHED_CLUSTER_PATH'] = self.check_dir('%s/published/' % os.environ.get('ODY_SEQ_ROOT', DEFAULT_SEQ_ROOT).rstrip('/'))
        self.data['REF_PATH'] = self.check_dir('%s/' % os.environ.get('ODY_REF', '%s/refs/10x/2019.05.19/cellranger' % DEFAULT_INFORMATICS_ROOT).rstrip('/'))
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

    def check_dir(self, path):
        if os.path.isdir(path):
            return path
        else:
            raise Exception('path does not exist: %s' % path)
