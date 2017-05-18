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
from subprocess import Popen,PIPE
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
    with readstream as stream:
        reads=grouper(stream,4)
        for read in reads:        
            prob_sampling=random.uniform(0,1)
            if prob_sampling<=samplingfreq and counter < max_sample_count:
                head,seq,placeholder,qual=[i.strip() for i in read]
        	counter+=1
                fastaout.write('>%s\n%s\n' % (head.split()[0],seq))
    
    fastaout.close()
    return fastaname
    

def dustmasker_runner(fasta,windowsize=20):             	
    dustcmd = 'dustmasker -in %s -out dustmasked_%s -window %s' % (fasta,fasta,windowsize) 
    dust_run = Popen(dustcmd,shell=True,stderr=PIPE,stdout=PIPE)
        dust_out,dust_err=dust_run.communicate() 
        if dust_run.returncode!=0:
            raise Exception('dustmasker run on %s failed: %s' % (fasta,dust_err)               

    return 'dustmasked_%s' % fasta       

def dust_summarize(fastqin,windowsize=20,samplingfreq=0.01,max_sample_count=100000):

    sampled_as_fasta = fastq_sampler(fastqin,samplingfreq=samplingfreq,max_sample_count=max_sample_count)
    dustmasked_fasta = dustmasker_runner(sampled_as_fasta,windowsize=windowsize)

    masked_proportions = []
    readstream = get_input_stream(fqin)
    counter=0
    with readstream as stream:
        reads=grouper(stream,4)
        for read in reads:
            head,seq,placeholder,qual=[i.strip() for i in read]
            
    masked_proportions.append(seq.count('N')/float(len(seq)))

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
    filestring = basename(fqin.replace('.fastq','').replace('.gz','')
    report = open('%s_entropy.txt' % filestring,'w')
    for key in dust_summary.keys()[:-1]:
        report.write('%s:\t%s\n' % (key,dust_summary_keys[key])
    report.close()
    
    vals_R_vector='vals=c(%s)' % ','.join([str(val) for val in dust_summary['vector']])
    Rcode = open('%s_dustmaskerhist.R' % filestring,'w')
    Rcode.write('%s\n' % vals_R_vector)
    outpdf = '%s.dustmasker.pdf' % filestring
    pdfcmd = 'pdf(file="%s.pdf")' % outpdf
    plotcmd = 'hist(vals,breaks=20,freq=FALSE,xlab=\"Low complexity proportion by read\")'
    closecmd = 'dev.off()'
    for cmd in [outpdf,pdfcmd,plotcmd,closecmd]:
        Rcode.write('%s\n' % plotcmd)
    
    makeplot = 'Rscript %s_dustmaskerhist.R' % filestring
    plot_run = Popen(makeplot,shell=True,stderr=PIPE,stdout=PIPE)
        plotrun_out,plotrun_err=plot_run.communicate() 
        if plot_run.returncode!=0:
            raise Exception('Error generating dustmasker histogram for %s: %s\n' % (fqin,plotrun_err))

    Rcode.close()



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="options for generating dustmasker metrics from fastq_file")
    parser.add_argument('-i','--fastqin',dest='fqin',type=str,help='fastq infile')
    parser.add_argument('-f','--subsample_freq',default=0.01,dest='subsamp',type=float,help='proportion of sequences to sample')
    parser.add_argument('-m','--max_samples',default=100000,dest=maxsamp,type=int,help='max number of reads to run dustmasker on')
    parser.add_argument('-w','--dust_window_size',default=20,dest=winsize,type=int,help='dustmasker windowsize')
    opts = parser.parse_args()

    write_dust_report(opts.fqin,windowsize=opts.winsize,samplingfreq=opts.subsamp,max_sample_count=opts.maxsamp)
    
