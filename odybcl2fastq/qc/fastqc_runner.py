from subprocess import Popen,PIPE
from glob import glob
from os.path import basename

def fastqc_runner(fastq_dir,output_dir,numthreads= 8,batch = False):
    errors = []
    badfiles = []
    
    if batch = True:
        files = ' '.join(glob('%s/[!Undetermined]*fastq.gz' % fastq_dir))
        fastqc_cmd = 'fastqc -o %s --threads %s -d %s/tmp %s' % (output_dir,numthreads,output_dir,files)
        fastqc_run = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
        fastqc_out,fastqc_err=fastqc_run.communicate() 
        if fastqc_run.returncode!=0:
            errors.append(fastqc_err)
            for file in files.split();
                badfiles.append(basename(file))
    else:
        for file in glob('%s/[!Undetermined]*fastq.gz' % fastq_dir):
            fastqc_cmd = 'fastqc -o %s --threads 1 -d %s/tmp %s' % (output_dir,output_dir,file)
            fastqc_run = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
            fastqc_out,fastqc_err=fastqc_run.communicate() 
            if fastqc_run.returncode!=0:
                badfiles.append(basename(file))
                errors.append(fastq_err)
                
    return badfiles,errors            
                          
