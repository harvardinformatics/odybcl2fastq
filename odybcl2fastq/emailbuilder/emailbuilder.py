import smtplib
from jinja2 import Environment, PackageLoader
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, make_msgid
from odybcl2fastq import config


def generateMessageId():
    '''
    Creates a unique message id for the email message
    '''
    return make_msgid()

def composeMessage(message, subject, summary_data, template, fromaddr, toemaillist, ccemaillist=[], bccemaillist=[]):
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
        html = get_html(summary_data, template)
        msg.attach(MIMEText(html.encode('utf-8'), 'html'))
    else:
        msg.attach(MIMEText(message.encode('utf-8')[-900000:], 'plain'))
    return msg

def buildmessage(message, subject, summary_data, template, fromaddr, toemaillist, ccemaillist=[], bccemaillist=[], server=None):
    if not server:
        server = config.EMAIL['smtp']
    msg = composeMessage(message, subject, summary_data, template, fromaddr, toemaillist, ccemaillist=[], bccemaillist=[])
    emails = toemaillist + ccemaillist + bccemaillist
    smtp = smtplib.SMTP(server)
    success = smtp.sendmail(fromaddr, emails, msg.as_string())
    smtp.close()
    return success

def get_html(summary_data, template):
    # create html message with jinja
    j2_env = Environment(
        loader=PackageLoader('odybcl2fastq', 'templates'),
        trim_blocks=True
    )
    html = j2_env.get_template(template).render(summary_data)
    return html
