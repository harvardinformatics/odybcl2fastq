# odybcl2fastq
A python packages wrapping Illumina's bcl2fastq software to demultiplex HiSeq, MiSeq, NextSeq, and NovaSeq runs, for use on the [FASRC cluster](https://www.rc.fas.harvard.edu/cluster/).
The end goal of this package is to automate demultiplexing of sequencing runs, and generate QC metrics and reports for end users that will aid in interpreting the quality of the resulting sequence data.

## Getting Started

1. Clone this git repository.
2. Download [bcl2fastq2 Conversion Software v2.20 Installer (Linux rpm)(https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html), convert the rpm to .deb using [alien](https://help.ubuntu.com/community/RPM/AlienHowto), and place in the root of the git working tree.
3. Build the Singularity image:

    singularity build ody_YYYYMMDD.sif ody_recipe.prod
4. Build separate images for cellranger and cellranger-atac (to be bind-mounted in the ody container at /opt/cellranger and /opt/cellranger-atac, respectively); e.g.:

```sh
mkdir -p data/opt
tar -C data/opt -xzf cellranger-atac-1.2.0.tar.gz
/usr/sbin/mksquashfs data cellranger-atac-1.2.0.squashfs -comp xz -Xbcj x86 -Xdict-size 1M -all-root
rm -rf data/opt && mkdir data/opt
tar -C data/opt -xzf cellranger-4.0.0.tar.gz
/usr/sbin/mksquashfs data cellranger-4.0.0.squashfs -b 1M -comp xz -Xbcj x86 -Xdict-size 1M -all-root
```

### Configuration

All configuration is set via environment variables and bind mounts (except for snakemake Slurm profile in odybcl2fastq/profiles/rc_slurm).


## Running Odybcl2fastq

### Multiple Runs
Use the script odybcl2fastq/process_snakemake_runs.py
This script will search the configured directory for runs which all the required
files and will execute a snakemake workflow for each run.


## Odybcl2fastq Logging

### Multiple Runs Logging
odybcl2fastq/process_snakemake_runs.py will log all the runs it queues to the
/log/odybcl2fastq10x.log.  This log also reports the sucess or failure of those runs.


### Single Run Log
Each run also gets it's own log file in /log/<run>.
This log will show the bcl2fastq cmd run as well as any output
from that job and post process jobs.

## Odybcl2fastq Alerting

### Multiple Run Alerting
An email is sent if a run fails, or an exception if encountered


### Single Run Alerting
An email is sent for any failure.  A warning is sent if outputdir space is close
to capacity.

For successfull runs an html summary of results is emailed to the address.


## Testing
Set di
There is a docker container, Dockerfile-test, that can be used to run tests in a CentOS 6 environment
with bcl2fastq installed (downloaded directly from Illumina).

To run the odybcl2fastq_tests.py, you'll need to map in the test illumina run to tests/test_run and
set the config file via ODYBCL2FASTQ_CONFIG_FILE to tests/test.config.json.  Also, a log file must be set.
