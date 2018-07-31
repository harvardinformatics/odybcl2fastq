from odybcl2fastq import constants as const
from odybcl2fastq import util
import os

CONFIG_FILE = os.environ.get('ODYBCL2FASTQ_CONFIG_FILE', os.path.join(const.APP_DIR, 'config.json'))


class Config(object):

    def __init__(self):
        config = util.load_json(CONFIG_FILE)
        for key, vars in config.items():
            setattr(self, key, vars)
