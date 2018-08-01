from odybcl2fastq.config import Config
import logging, os
from odybcl2fastq import constants as const


__version__ = '1.0.3'

config = Config()

LOGLEVELSTR = os.environ.get('ODYBCL2FASTQ_LOG_LEVEL', 'INFO')

logging.basicConfig(
    level=logging.getLevelName(LOGLEVELSTR),
    format='%(asctime)s %(message)s'
)

# Setup loggers
logger  = logging.getLogger('odybcl2fastq')
logfile = os.environ.get('ODYBCL2FASTQ_LOG_FILE', os.path.join(const.ROOT_DIR, 'odybcl2fastq.log'))
handler = logging.FileHandler(logfile)
handler.setLevel(logging.getLevelName(LOGLEVELSTR))
logger.addHandler(handler)

logger  = logging.getLogger('centrifuge')
logfile = os.environ.get('CENTRIFUGE_LOG_FILE', os.path.join(const.ROOT_DIR, 'centrifuge.log'))
handler = logging.FileHandler(logfile)
handler.setLevel(logging.getLevelName(os.environ.get('CENTRIFUGE_LOG_LEVEL', LOGLEVELSTR)))
logger.addHandler(handler)

logger  = logging.getLogger('storage_mgt')
logfile = os.environ.get('STORAGEMGT_LOG_FILE', os.path.join(const.ROOT_DIR, 'storage_mgmt.log'))
handler = logging.FileHandler(logfile)
handler.setLevel(logging.getLevelName(os.environ.get('STORAGEMGT_LOG_LEVEL', LOGLEVELSTR)))
logger.addHandler(handler)

logger  = logging.getLogger('db')
logfile = os.environ.get('DB_LOG_FILE', os.path.join(const.ROOT_DIR, 'db.log'))
handler = logging.FileHandler(logfile)
handler.setLevel(logging.getLevelName(os.environ.get('DB_LOG_LEVEL', LOGLEVELSTR)))
logger.addHandler(handler)


class UserException(Exception):
    '''
    Exceptions for known user errors
    '''

    def __init__(self, message):
        super(UserException, self).__init__(message)
        self.user_msg = message
