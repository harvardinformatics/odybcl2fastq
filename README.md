# odybcl2fastq
A python packages wrapping Illumina's bcl2fastq software to demultiplex HiSeq and NextSeq runs, for use on the Odyssey cluster. The end goal of this package is to automate demultiplexing of sequencing runs, and generate QC metrics and reports for end users that will aid in interpreting the quality of the resulting sequence data. Odybcl2fastq is being developed using the Anaconda 1.9.2 distribution of python (v. 2.7.11).
## Getting Started

Clone this git repository.  Set up an anaconda environment with python 2.7.11. Install the dependencies.

### Dependencies

* jinja2: conda install jinja2
* bcl2fastq: this can be obtained from illumina, bcl2fastq2 version 2.2
* fastqc: fastqc version 0.11.5

### Config

Set the variables in config.json.  See the example config.json.example, insert
values for all the variables and then save the file as config.json.


## Running Odybcl2fastq

### Single Run
Use the script odybcl2fastq/run.py
Required arguments:
* --runinfoxml: path to the run info file
* --sample-sheet: path to the sample sheet file
* --runfolder: path to the run folder
* --output-dir: path to the output

Optional arguments:
* --test: prints out the bcl2fastq cmd but does not run it, then exits
* --no-demultiplex: skips the demultiplexing part of the script
* --no-post-process: skips updating the lims db and running fastqc
* --no-file-copy: skips copying from output dir to final dir
Many bcl2fastq parameters are options for a full list please see parameter defs
in odybcl2fastq/run.py file


### Multiple Runs
Use the script odybcl2fastq/process_runs.py
This script will search the configured directory for runs which all the required
files and will queue off a pool of odybcl2fastq/run.py calls for each run.

Environment variable, 'ODYBCL2FASTQ_PROC_NUM', will determine the number of
parellel runs to processs.


## Odybcl2fastq Logging

### Multiple Runs Logging
odybcl2fastq/process_runs.py will log all the runs it queues to the
odybcl2fastq.log.  This log also reports the sucess or failure of those runs.


### Single Run Log
Each run also gets it's own log file, the location of these is configured in
config.json.  This log will show the bcl2fastq cmd run as well as any output
from that job and post process jobs.

## Odybcl2fastq Alerting

### Multiple Run Alerting
An email is sent if a run fails, or an exception if encountered


### Single Run Alerting
An email is sent for any failure.  A warning is sent if outputdir space is close
to capacity.

For successfull runs an html summary of results is emailed.


## Testing
There is a docker container, Dockerfile-test, that can be used to run tests in a CentOS 6 environment
with bcl2fastq installed (downloaded directly from Illumina).

To run the odybcl2fastq_tests.py, you'll need to map in the test illumina run to tests/test_run and
set the config file via ODYBCL2FASTQ_CONFIG_FILE to tests/test.config.json.  Also, a log file must be set.

e.g.

    docker build -t odybcl2fastq -f Dockerfile-test .
    docker run -it -e ODYBCL2FASTQ_CONFIG_FILE=/app/tests/test.config.json -v $HOME/test/odybcl2fastq/test_run:/app/tests/test_run -e ODYBCL2FASTQ_LOG_LEVEL=DEBUG -e ODYBCL2FASTQ_LOG_FILE=/app/logs/odyfile.log odybcl2fastq nosetests -vv tests/odybcl2fastq_tests.py
