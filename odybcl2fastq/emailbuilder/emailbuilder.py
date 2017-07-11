import smtplib
from jinja2 import Environment, FileSystemLoader
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE,make_msgid
import constants as const

def generateMessageId():
    '''
    Creates a unique message id for the email message
    '''
    return make_msgid()


def buildmessage(message, subject, lane_data, sample_data, fromaddr,toemaillist, ccemaillist=[], bccemaillist=[], server='smtp.fas.harvard.edu'):
    msg = MIMEMultipart()
    msg['Message-ID'] = generateMessageId()
    msg['From'] = fromaddr
    msg['To'] = COMMASPACE.join(toemaillist)
    msg['Subject'] = subject
    if len(ccemaillist) > 0:
        msg['Cc'] = COMMASPACE.join(ccemaillist)
    if len(bccemaillist) > 0:
        msg['Bcc'] = COMMASPACE.join(bccemaillist)
    html = get_html(message, subject, lane_data, sample_data)
    msg.attach(MIMEText(html.encode('utf-8'),'html'))
    emails = toemaillist + ccemaillist + bccemaillist
    smtp = smtplib.SMTP(server)
    smtp.sendmail(fromaddr,emails,msg.as_string())
    smtp.close()

def get_html(run, message, lane_data, sample_data):
    # create html message with jinja
    j2_env = Environment(loader=FileSystemLoader(const.TEMPLATE_DIR),
            trim_blocks = True)
    # context to put in email template
    context = {
            'run': run,
            'subject': run,
            'analyses': 1,
            'sample_no': len(sample_data),
            'reads': ['363','058','633'],
            'lane_data': lane_data,
            'sample_data': sample_data,
            'letter': '',
            # TODO: store in config or something
            'run_folder': 'http://software.rc.fas.harvard.edu/ngsdata/'

    }
    html = j2_env.get_template('summary.html').render(context)
    return html


