import os
import gzip
from glob import glob
from odybcl2fastq import config
import json
import logging

def compare_fastq(test_dir, instrument, run):
    """
    parse summary from Stats.json
    """
    errors = []
    control_dir = config.CONTROL_DIR
    if not os.path.exists(control_dir):
        raise UserException('control dir does not exist: %s' % control_dir)
    file_match = '/*/Fastq/*.fastq.gz'
    control_dir += run
    control_fastq = glob(control_dir + file_match)
    files_checked = []
    for f in control_fastq:
        filename = f.split('/')[-1].split('.')[0].replace('-', '_')
        test_match = '/[!QC]*/' + filename + '*'
        test_fastq = glob(test_dir + test_match)
        test_match = '/' + filename + '*'
        test_path = test_dir + test_match
        t_glob = glob(test_path)
        test_fastq.extend(t_glob)
        if not test_fastq:
            errors.append('skipping this file since cant find test %s' % test_match)
        elif len(test_fastq) != 1:
            errors.append('more than one file match %s' % test_match)
        else:
            files_checked.append(filename)
            test_file = test_fastq[0]
            c = gzip.open(f)
            t = gzip.open(test_file)
            for i, c_line in enumerate(c):
                t_line = t.readline()
                if i % 4 == 0:
                    c_header = c_line.split(':')[2:7]
                    t_header = t_line.split(':')[2:7]
                    if c_header != t_header:
                        errors.append('%s differences in header %s vs %s' %
                                (filename, json.dumps(c_header), json.dumps(t_header)))
                elif i % 4 == 1:
                    if instrument == 'hiseq':
                        c_line = reverse_complement(c_line)
                    if not lines_match(c_line, t_line):
                        errors.append('%s line %i differences in seq control (%i): %s vs test (%i): %s' %
                                (filename, i, len(c_line), c_line, len(t_line), t_line))
                elif i % 4 == 3:
                    if instrument == 'hiseq':
                        c_line = reverse(c_line)
                    if not lines_match(c_line, t_line):
                        errors.append('%s differences in score %s vs %s' %
                                (filename, c_line, t_line))
                if i == 400:
                    break
        if len(files_checked) == 3:
            break
    print(files_checked)
    return errors

def lines_match(c_line, t_line):
    c_line = c_line.replace('/n', '').strip()
    t_line = t_line.replace('/n', '').strip()
    t_len = len(t_line)
    for k, char in enumerate(c_line):
        if k < t_len:
            if char != t_line[k]:
                return False
        else:
            break
    return True

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])

def reverse(seq):
    return ''.join([base for base in seq[::-1]])

