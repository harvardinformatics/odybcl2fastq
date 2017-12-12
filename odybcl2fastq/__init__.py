from odybcl2fastq.config import Config

config = Config()

class UserException(Exception):
    '''
    Exceptions for known user errors
    '''

    def __init__(self,message):
        super(UserException,self).__init__(message)
        self.user_msg = message
