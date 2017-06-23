'''
randomly samples a fastq file, converts reads to fasta,
then runs dustmasker on it to generate a low-complexity masked
fasta from which histogram and stats on complexity are generated

requires loading ncbi tools module:
   module load gcc/4.8.2-fasrc01 ncbi-tools++/12.0.0-fasrc01
'''

import argparse
import gzip
from itertools import izip,izip_longest
from numpy import mean,median,log2,random,std,histogram
from scipy.stats import sem
from collections import defaultdict,OrderedDict
from os.path import basename
from subprocess import Popen,PIPE,call
from odybcl2fastq import UserException


def get_input_stream(fastqfile):
    if fastqfile[-2:]=='gz':
        fastqhandle=gzip.open(fastqfile,'rb')
    else:
        fastqhandle=open(fastqfile,'r')
    
    return fastqhandle    


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def fastq_sampler(fqin,samplingfreq=0.01,max_sample_count=100000):
    fastaname = '%s.fasta' % basename(fqin).replace('.fastq','').replace('.gz','')
    fastaout=open(fastaname,'w')
    readstream = get_input_stream(fqin)
    counter=0
    read_lengths={}
    with readstream as stream:
        reads=grouper(stream,4)
        for read in reads:        
            prob_sampling=random.uniform(0,1)
            if prob_sampling<=samplingfreq:
                head,seq,placeholder,qual=[i.strip() for i in read]
                head = head.split()[0]                
        	counter+=1
                fastaout.write('>%s\n%s\n' % (head.split()[0],seq))
                read_lengths[head[1:]] = len(seq)
            if counter == max_sample_count:
                break
    
    fastaout.close()
    return fastaname,read_lengths
    

def dustmasker_runner(fasta,windowsize=20):             	
    dustcmd = 'dustmasker -in %s -out dustmasked_%s.intervals -window %s -outfmt acclist' % (fasta,fasta,windowsize) 
    dust_run = Popen(dustcmd,shell=True,stderr=PIPE,stdout=PIPE)
    dust_out,dust_err=dust_run.communicate() 
    if dust_run.returncode!=0:
        raise Exception('dustmasker run on %s failed: %s' % (fasta,dust_err))               
    else:
        return 'dustmasked_%s.intervals' % fasta       

def dust_summarize(fastqin,windowsize=20,samplingfreq=0.01,max_sample_count=100000):
    masked_proportions = []
    sampled_as_fasta,read_lengths = fastq_sampler(fastqin,samplingfreq=samplingfreq,max_sample_count=max_sample_count)
    fasta_open=open(sampled_as_fasta,'r')
    print 'read lengths', read_lengths   
    dustmasked_intervals = dustmasker_runner(sampled_as_fasta,windowsize=windowsize)
    read_mask_dict={}
    interval_parse=open(dustmasked_intervals,'r')
    for line in interval_parse:
        linelist=line.replace('>@','').strip().split('\t')
    
        read_mask_dict[linelist[0]] = (int(linelist[2]) - int(linelist[1]) + 1)/float(read_lengths[linelist[0]])
    for read in read_lengths.keys():
        if read not in read_mask_dict:
            masked_proportions.append(0)
        else:
            masked_proportions.append(read_mask_dict[read])
        

    dust_summary = OrderedDict()
    dust_summary['nsamples'] = len(masked_proportions)
    dust_summary['mean'] = mean(masked_proportions)
    dust_summary['median'] = median(masked_proportions)
    dust_summary['std'] = std(masked_proportions)
    dust_summary['sem'] = sem(masked_proportions)
    dust_summary['vector'] = masked_proportions
        
    return dust_summary

def write_dust_report(fqin,windowsize=20,samplingfreq=0.01,max_sample_count=100000):
    dust_summary=dust_summarize(fqin,windowsize=windowsize,samplingfreq=samplingfreq,max_sample_count=max_sample_count)
    filestring = basename(fqin.replace('.fastq','').replace('.gz',''))
    report = open('%s_dustmasker.txt' % filestring,'w')
    for key in dust_summary.keys()[:-1]:
        report.write('%s:\t%s\n' % (key,dust_summary[key]))
    report.close()
    
    vals_R_vector='vals=c(%s)' % ','.join([str(val) for val in dust_summary['vector']])
    Rcode = open('%s_dustmaskerhist.R' % filestring,'w')
    Rcode.write('%s\n' % vals_R_vector)
    outpdf = '%s.dustmasker.pdf' % filestring
    pdfcmd = 'pdf(file="%s")' % outpdf
    plotcmd = 'hist(vals,breaks=20,freq=FALSE,main="",xlab=\"Low complexity proportion by read\")'
    closecmd = 'dev.off()'
    for cmd in [pdfcmd,plotcmd,closecmd]:
        Rcode.write('%s\n' % cmd)
    Rcode.close()

    makeplot = 'Rscript %s_dustmaskerhist.R' % filestring
    call("%s" % makeplot,shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="options for generating dustmasker metrics from fastq_file")
    parser.add_argument('-i','--fastqin',dest='fqin',type=str,help='fastq infile')
    parser.add_argument('-f','--subsample_freq',default=0.01,dest='subsamp',type=float,help='proportion of sequences to sample')
    parser.add_argument('-m','--max_samples',default=100000,dest='maxsamp',type=int,help='max number of reads to run dustmasker on')
    parser.add_argument('-w','--dust_window_size',default=20,dest='winsize',type=int,help='dustmasker windowsize')
    opts = parser.parse_args()

    write_dust_report(opts.fqin,windowsize=opts.winsize,samplingfreq=opts.subsamp,max_sample_count=opts.maxsamp)
    
