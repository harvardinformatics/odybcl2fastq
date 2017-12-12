from odybcl2fastq import constants as const
from odybcl2fastq import util

CONFIG_FILE = 'config.json'

class Config(object):

    def __init__(self):
        config = util.load_json(const.APP_DIR + CONFIG_FILE)
        for key, vars in config.items():
            setattr(self, key, vars)
