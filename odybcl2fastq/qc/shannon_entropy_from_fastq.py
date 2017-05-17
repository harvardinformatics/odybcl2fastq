import argparse
import gzip
from itertools import izip,izip_longest
from numpy import mean,median,log2,random,std,histogram
from scipy.stats import sem
from collections import defaultdict,OrderedDict
from sets import Set
from os.path import basename
import matplotlib.pyplot as plt
from subprocess import Popen,PIPE

def get_input_stream(fastqfile):
    if fastqfile[-2:]=='gz':
        fastqhandle=gzip.open(fastqfile,'rb')
    else:
        fastqhandle=open(fastqfile,'r')
    
    return fastqhandle    


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
    

def stripNs(rawsequence):
    cleanseq=rawsequence.replace('N','')
    return cleanseq
    
    
def calc_base_freqs(stripped_sequence):
    basecounts=defaultdict(int)
    for base in Set(stripped_sequence):
        basecounts[base]+=stripped_sequence.count(base)
    basefreqs={}
    for base in basecounts.keys():
       basefreqs[base]=basecounts[base]/float(sum(basecounts.values()))
       
    return basefreqs    

def shannon_entropy(base_freq_dict):
    sentropy=-1*sum([log2(freq) for freq in base_freq_dict.values()])
    return sentropy
    

def entropy_runner(fqin,samplingfreq,max_sample_count = 100000):
    readstream = get_input_stream(fqin)
    entropy_values = []
    
    with readstream as stream:
        reads=grouper(stream,4)
        while len(entropy_values) <= max_sample_count:
            for read in reads:        
                prob_sampling=random.uniform(0,1)
                if prob_sampling<=samplingfreq and len(enotropy_values)<=max_sample_count:
                    head,seq,placeholder,qual=[i.strip() for i in read]
                    nstrip_read=stripNs(seq)
                    nucleotide_freqs=calc_base_freqs(nstrip_read)
                    read_entropy=shannon_entropy(nucleotide_freqs)
                    entropy_values.append(read_entropy)
        	    	

        entropy_summary=OrderedDict()
        entropy_summary['nsamples'] = len(entropy_values)
        entropy_summary['mean'] = mean(entropy_values)
        entropy_summary['median'] = median(entropy_values)
        entropy_summary['std'] = std(entropy_values)
        entropy_summary['sem'] = sem(entropy_values)
        entropy_summary['vector'] = entropy_values   
        

    return entropy_summary

def write_entropy_report(fqin,samplingfreq,**kwargs):
    entropy_summary=entropy_runner(fqin,samplingfreq,**kwargs)
    filestring = basename(fqin.replace('.fastq','').replace('.gz','')
    report = open('%s_entropy.txt' % filestring,'w')
    for key in entropy_summary.keys()[:-1]:
        report.write('%s:\t%s\n' % (key,entropy_summary_keys[key])
    report.close()
    
    vals_R_vector='vals=c(%s)' % ','.join([str(val) for val in entropy_summary['vector']])
    Rcode = open('%s_entropyhist.R' % filestring,'w')
    Rcode.write('%s\n' % vals_R_vector)
    outpdf = '%s.entropy.pdf' % filestring
    pdfcmd = 'pdf(file="%s.pdf")' % outpdf
    plotcmd = 'hist(vals,breaks=20,freq=FALSE,xlab=\"Entropy\")'
    closecmd = 'dev.off()'
    for cmd in [outpdf,pdfcmd,plotcmd,closecmd]:
        Rcode.write('%s\n' % plotcmd)
    
    makeplot = 'Rscript %s_entropyhist.R' % filestring
    plot_run = Popen(makeplot,shell=True,stderr=PIPE,stdout=PIPE)
        plotrun_out,plotrun_err=plot_run.communicate() 
        if plot_run.returncode!=0:
            print('Error generating entropy histogram for %s: %s\n' % (fqin,plotrun_err))

    Rcode.close()
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="options for calculating entropy from a fastq_file")
    parser.add_argument('-i','--fastqin',dest='fqin',type=str,help='fastq infile')
    parser.add_argument('-f','--subsample_freq',default=0.001,dest='subsamp',type=float,help='proportion of sequences to print entropy for')
    parser.add_argument('-m','--max_samples',default=100000,dest=maxsamp,type=int,help='max number of reads to calculate entropy')
    opts = parser.parse_args()
    
    write_entropy_report(opts.fqin,opts.subsamp,opts.maxsamp)
