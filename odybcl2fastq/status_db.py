import os
import json
import time
import MySQLdb
from odybcl2fastq import config

class StatusDB(object):
    def __init__(self):
        self.db = MySQLdb.connect(
                host = config.STATUSDB['host'],
                user = config.STATUSDB['user'],
                passwd = config.STATUSDB['password'],
                db = config.STATUSDB['database']
        )

    def minilims_select(self, thing = None, name = None, property = None, value
            = None):
        # get dict of func args
        args = locals()
        # valid columns
        cols = ['thing', 'name', 'property', 'value']
        if not thing and not name and not property and not value:
            raise UserException('no criteria for minilims query')
        sql = 'select * from semantic_data where '
        wheres = []
        # add args to wheres
        for k, v in args.items():
            if k in cols and v:
                wheres.append("%s = '%s'" % (k, v))
        sql += ' and '.join(wheres)
        cur = self.db.cursor()
        rows = cur.execute(sql)
        rows = cur.fetchall()
        cur.close()
        return rows

    def minilims_get_new_name(self, table, prefix):
        sql = """select name from semantic_data where thing = '%s'
            and property = 'name' order by name desc limit 1""" % table
        cur = self.db.cursor()
        cur.execute(sql)
        prev_name = cur.fetchone()[0]
        cur.close()
        prev_num = int(prev_name.split(prefix)[1])
        next_num = str(prev_num + 1)
        num_char = len(next_num)
        for i in range(num_char, 5):
            next_num = '0' + next_num
        return prefix + next_num

    def minilims_insert(self, data, table, name):
        cur = self.db.cursor()
        for k, v in data.items():
            sql = """insert into semantic_data (name, thing, property, value) values ('%s',
            '%s', '%s', '%s')""" % (name, table, k, v)
            cur.execute(sql)
        self.db.commit()
        cur.close()

    def insert_analysis(self, run, subs_str):
        analysis_table = 'Illumina_BclConversion_Analysis'
        rows = self.minilims_select(thing=analysis_table, property='Illumina_Run', value=run)
        row_cnt = len(rows)
        now = time.time()
        if row_cnt == 0:
            name = self.minilims_get_new_name(analysis_table, 'ILL')
            data = {
                'Data_Directory': config.MOUNT_DIR + run,
                'Date_Created': now,
                'Date_Modified': now,
                'Illumina_BclConversion_Analysis': name,
                'Illumina_Run': run,
                'Name': name,
                'Status': 'COMPLETE',
                'Submission': subs_str,
                'Web_Link': config.FASTQ_URL + run
            }
            self.minilims_insert(data, 'Illumina_Bcl_Conversion_Analysis', name)
        else:
            name = rows[0]['name']
            if row_cnt > 1: # dups exists
                logging.warning("duplicate analysis: %s\n", name)
        self.analysis_name = name

    def link_run_and_subs(self, run, subs):
        for sub in subs:
            self.minilims_insert({'Illumina_run': run}, 'Submission',  sub)
            self.minilims_insert({'Submission': sub}, 'Illumina_Run', run)


