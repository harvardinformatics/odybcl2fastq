import argparse
from subprocess import Popen,PIPE
from glob import glob
import os
from subprocess import call
import gzip
import logging

def fastqc_runner(output_dir,numthreads = 1,batch = False):
    errors = []
    badfiles = []
    # get fastq files except undetermined
    files = ' '.join(glob('%s/*/[!Undetermined]*fastq*.gz' % output_dir))
    print 'files are', files
    # create qc dir
    qc_dir = output_dir + '/QC'
    if not os.path.exists(qc_dir):
        call('mkdir %s' % (qc_dir) ,shell=True)
    if batch == True:
        fastqc_cmd = 'fastqc -o %s --threads %s -b %s' % (qc_dir,numthreads,files)
        fastqc_run = Popen(fastqc_cmd,shell=True,stderr=PIPE,stdout=PIPE)
        fastqc_out,fastqc_err=fastqc_run.communicate()

        if fastqc_run.returncode!=0:
            errors.append(fastqc_err)
            for file in files.split():
                badfiles.append(os.pathbasename(file))
            out.append(fastqc_out)
    else:
        for file in files.split():
            print 'file is', os.path.basename(file)
            if gzNotEmpty(file):
                print('gz is not empty' + file)
                fastqc_cmd = 'fastqc -o %s --noextract --threads 1 %s' % (qc_dir,file)
                logging.info("FASTQC: " + fastqc_cmd)
                fastqc_run = Popen(fastqc_cmd,shell=True,stderr=PIPE,stdout=PIPE)
                fastqc_out,fastqc_err=fastqc_run.communicate()
                if fastqc_run.returncode!=0:
                    logging.info("FASTQC failed for: " + file)
                else:
                    logging.info("FASTQC complete successfully for: " + file)
                    print(fastqc_out)
                    print(fastqc_err)
                badfiles.append(os.path.basename(file))
                errors.append(fastqc_err)
                out.append(fastqc_out)
    return badfiles, errors, out

def gzNotEmpty(file):
    fh = gzip.open(file, 'rb')
    data = fh.read(100)
    fh.close()
    return bool(data)

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
