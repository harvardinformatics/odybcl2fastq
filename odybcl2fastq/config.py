import json
import os
import subprocess

class Config(object):

    def __init__(self):
        self.data = {}
        self.data['EMAIL_ADMIN'] = json.loads(os.environ['ODY_EMAIL_ADMIN'])
        self.data['EMAIL_FROM'] = os.environ['ODY_EMAIL_FROM']
        self.data['EMAIL_SMTP'] = os.environ['ODY_EMAIL_SMTP']
        self.data['EMAIL_TO'] = json.loads(os.environ['ODY_EMAIL_TO'])
        self.data['FASTQ_URL'] = os.environ.get('ODY_FASTQ_URL', 'https://software.rc.fas.harvard.edu/ngsdata/')
        self.data['GLOBUS_URL'] = os.environ['ODY_GLOBUS_URL']
        self.data['PUBLISHED_CLUSTER_PATH'] = os.environ['ODY_PUBLISHED_CLUSTER_PATH']
        self.data['ANALYSIS_DIR'] = os.environ['ODY_ANALYSIS_DIR']
        # pre-flight checks to ensure existence and accessibility of required directories
        self.check_dir('/sequencing/source')
        self.check_dir('/sequencing/analysis', check_is_writable=True)
        self.check_dir('/sequencing/published', check_is_writable=True)
        self.check_dir('/ref')
        self.data['TEST'] = os.environ.get('ODY_TEST', 'FALSE') == 'TRUE'
        # default to empty for db connection variables for now since they are
        # not used in when TEST is true
        self.data['STATUS_DB_HOST'] = os.environ.get('ODY_STATUS_DB_HOST', '')
        self.data['STATUS_DB_NAME'] = os.environ.get('ODY_STATUS_DB_NAME', '')
        self.data['STATUS_DB_USER'] = os.environ.get('ODY_STATUS_DB_USER', '')
        self.data['STATUS_DB_PASSWORD'] = os.environ.get('ODY_STATUS_DB_PASSWORD', '')
        self.data['BAUER_API'] = os.environ.get('ODY_BAUER_API', '')
        self.data['BAUER_TOKEN'] = os.environ.get('ODY_BAUER_TOKEN', '')

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

    def check_dir(self, path, check_is_writable=False):
        if os.path.isdir(path):
            if not os.access(path, os.R_OK | os.X_OK):
                raise Exception('path is not readable: %s' % path)
            if check_is_writable and not os.access(path, os.W_OK):
                raise Exception('path is not writable: %s' % path)
            else:
                return path
        else:
            raise Exception('path does not exist: %s' % path)
