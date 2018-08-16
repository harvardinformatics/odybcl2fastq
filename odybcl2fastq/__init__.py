from odybcl2fastq.config import Config
import logging, os
from odybcl2fastq import constants as const


__version__ = '1.0.3'

config = Config()

LOGLEVELSTR = os.environ.get('ODYBCL2FASTQ_LOG_LEVEL', 'INFO')

logging.basicConfig(
    level=logging.getLevelName(LOGLEVELSTR),
    format='%(message)s'
)

# Setup loggers
logger  = logging.getLogger('odybcl2fastq')
logfilename = os.environ.get('ODYBCL2FASTQ_LOG_FILE', 'odybcl2fastq.log')
if not logfilename.startswith('/'):
    logfilename = os.path.join(config.LOG_DIR, logfilename)
handler = logging.FileHandler(logfilename)
handler.setLevel(logging.getLevelName(LOGLEVELSTR))
logger.addHandler(handler)

logger  = logging.getLogger('centrifuge')
logfile = os.environ.get('CENTRIFUGE_LOG_FILE', os.path.join(const.ROOT_DIR, 'centrifuge.log'))
handler = logging.FileHandler(logfile)
handler.setLevel(logging.getLevelName(os.environ.get('CENTRIFUGE_LOG_LEVEL', LOGLEVELSTR)))
logger.addHandler(handler)


def initLogger(name):
    '''
    Return a logger for the given name.  Uses the name to get log file and log level from
    the environment.
    If log file is not found in the env, stderr logging is used.
    If the log file is not an absolute path, the config.LOG_DIR is added.
    Formatter for log files is '%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p'
    '''
    logfileenv = '_'.join([name, 'log', 'file']).upper()
    loglevelenv = '_'.join([name, 'log', 'level']).upper()
    print "%s %s" % (logfileenv, loglevelenv)

    loglevel = logging.getLevelName(os.environ.get(loglevelenv, 'INFO'))
    logfilename = os.environ.get(logfileenv)

    logger = logging.getLogger(name)
    if logfilename:
        if not logfilename.startswith('/'):
            logfilename = os.path.join(config.LOG_DIR, logfilename)
        handler = logging.FileHandler(logfilename)
        handler.setLevel(loglevel)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        logger.setLevel(loglevel)

    return logger


class UserException(Exception):
    '''
    Exceptions for known user errors
    '''

    def __init__(self, message):
        super(UserException, self).__init__(message)
        self.user_msg = message
