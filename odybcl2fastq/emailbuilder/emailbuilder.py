import smtplib
from jinja2 import Environment, FileSystemLoader
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE,make_msgid
import odybcl2fastq.constants as const

def generateMessageId():
    '''
    Creates a unique message id for the email message
    '''
    return make_msgid()


def buildmessage(message, subject, summary_data, fromaddr,toemaillist, ccemaillist=[], bccemaillist=[], server='rcsmtp.rc.fas.harvard.edu'):
    msg = MIMEMultipart()
    msg['Message-ID'] = generateMessageId()
    msg['From'] = fromaddr
    msg['To'] = COMMASPACE.join(toemaillist)
    msg['Subject'] = subject
    if len(ccemaillist) > 0:
        msg['Cc'] = COMMASPACE.join(ccemaillist)
    if len(bccemaillist) > 0:
        msg['Bcc'] = COMMASPACE.join(bccemaillist)
    if summary_data:
        html = get_html(summary_data)
        msg.attach(MIMEText(html.encode('utf-8'),'html'))
    else:
        msg.attach(MIMEText(message.encode('utf-8'),'plain'))
    emails = toemaillist + ccemaillist + bccemaillist
    smtp = smtplib.SMTP(server)
    smtp.sendmail(fromaddr,emails,msg.as_string())
    smtp.close()

def get_html(summary_data):
    # create html message with jinja
    j2_env = Environment(loader=FileSystemLoader(const.TEMPLATE_DIR),
            trim_blocks = True)
    html = j2_env.get_template('summary.html').render(summary_data)
    return html


