import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE,make_msgid

def generateMessageId():
    '''
    Creates a unique message id for the email message
    '''
    return make_msgid()


def buildmessage(message,fromaddr,toemaillist,subject,ccemaillist=[],bccemaillist=[],server='rcsmtp.rc.fas.harvard.edu'):
    msg = MIMEMultipart()
    msg['Message-ID'] = generateMessageId()
    msg['From'] = fromaddr
    msg['To'] = COMMASPACE.join(toemaillist)
    msg['Subject'] = subject
    if len(ccemaillist) > 0:
        msg['Cc'] = COMMASPACE.join(ccemaillist)
    if len(bccemaillist) > 0:
        msg['Bcc'] = COMMASPACE.join(bccemaillist)

    msg.attach(MIMEText(message.encode('utf-8'),'html'))
    emails = toemaillist + ccemaillist + bccemaillist
    smtp = smtplib.SMTP(server)
    smtp.sendmail(fromaddr,emails,msg.as_string())
    smtp.close()    



         
