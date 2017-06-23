import argparse
from subprocess import Popen,PIPE
from glob import glob
from os.path import basename
from subprocess import call

def fastqc_runner(fastq_dir,output_dir,numthreads = 1,batch = False):
    errors = []
    badfiles = []
    files = ' '.join(glob('%s/[!Undetermined]*fastq*' % fastq_dir)) 
    print 'files are', files
    if batch == True:
        fastqc_cmd = 'fastqc -o %s --threads %s -d %s/tmp %s' % (output_dir,numthreads,output_dir,files)
        fastqc_run = Popen(fastqc_cmd,shell=True,stderr=PIPE,stdout=PIPE)
        fastqc_out,fastqc_err=fastqc_run.communicate() 
        if fastqc_run.returncode!=0:
            errors.append(fastqc_err)
            for file in files.split():
                badfiles.append(basename(file))
    else:
        for file in files.split():
            print 'file is', basename(file)
            call('mkdir %s/%s_tmp' % (output_dir,basename(file)),shell=True)
            fastqc_cmd = 'fastqc -o %s --threads 1 -d %s/%s_tmp %s' % (output_dir,output_dir,basename(file),file)
            fastqc_run = Popen(fastqc_cmd,shell=True,stderr=PIPE,stdout=PIPE)
            fastqc_out,fastqc_err=fastqc_run.communicate() 
            if fastqc_run.returncode!=0:
                badfiles.append(basename(file))
                errors.append(fastqc_err)
                
    return badfiles,errors            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="options for running fastqc")
    parser.add_argument('-o','--output_dir',dest='outdir',type=str,help='where to write fastqc output')
    parser.add_argument('-i','--fastq_dir',dest='indir',type=str,help='directory where fastq files are locate')
    parser.add_argument('-t','--threads',dest='nthreads',type=int,default=1,help="num files to process simultaneously")
    parser.add_argument('-b','--batch',dest='batch',type=str,default=False,help="whether to run as a batch")
    opts = parser.parse_args()

    badfiles,errors = fastqc_runner(opts.indir,opts.outdir,numthreads=opts.nthreads,batch=opts.batch)
    if len(errors) != 0:
        for i in range(len(errors)):
            print '%s\n%s\n' % (badfiles[i],errors[i])                     
