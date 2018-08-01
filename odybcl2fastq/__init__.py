from odybcl2fastq.config import Config

__version__ = '1.0.3'

config = Config()


class UserException(Exception):
    '''
    Exceptions for known user errors
    '''

    def __init__(self, message):
        super(UserException, self).__init__(message)
        self.user_msg = message
