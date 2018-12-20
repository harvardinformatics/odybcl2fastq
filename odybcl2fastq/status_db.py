import os
import json
import time
import MySQLdb
import logging
from odybcl2fastq import config

class StatusDB(object):
    def __init__(self):
        self.db = MySQLdb.connect(
                host = config.STATUS_DB['host'],
                user = config.STATUS_DB['user'],
                passwd = config.STATUS_DB['password'],
                db = config.STATUS_DB['database']
        )
        pass

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

    def minilims_insert(self, data, table, name, allow_dup = True):
        cur = self.db.cursor()
        for k, v in data.items():
            if allow_dup:
                dup = False
            else:
                dup = self.minilims_select(thing=table, name=name, property=k, value=v)
            if not dup:
                sql = """insert into semantic_data (name, thing, property, value) values ('%s',
                '%s', '%s', '%s')""" % (name, table, k, v)
                cur.execute(sql)
        self.db.commit()
        cur.close()

    def insert_analysis(self, run, subs_str):
        analysis_table = 'Illumina_BclConversion_Analysis'
        rows = self.minilims_select(thing=analysis_table, property='Illumina_Run', value=run)
        if len(rows) > 0:
            name = rows[0][0]
            logging.info("This run already has an analysis, %s this must be a rerun and an additional entry for this run will be created\n", name)
        now = time.time()
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
        self.minilims_insert(data, analysis_table, name)
        self.analysis_name = name
        return self.analysis_name

    def link_run_and_subs(self, run, subs):
        for sub in subs:
            # check if this is a rerun and this run has already been linked to
            # the submission
            existing_runs = self.minilims_select(thing='Submission', name=sub, property='Illumina_Run')
            exists = False
            for exist in existing_runs:
                if exist[3] in run:
                    exists = True
            allow_dup = False
            if not exists: # linking twice results in double billing
                self.minilims_insert({'Illumina_Run': run}, 'Submission',  sub,
                        allow_dup)
            # it is ok to link the submission to the any run it applys to
            self.minilims_insert({'Submission': sub}, 'Illumina_Run', run,
                    allow_dup)


