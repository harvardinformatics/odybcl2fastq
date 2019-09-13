#
# Docker container for dev odybcl2fastq
#
FROM centos:7.4.1708

RUN yum -y update && yum install -y python wget mysql-devel gcc make python-devel unzip libffi-devel libssl-devel pyOpenSSL
RUN yum -y install epel-release
RUN yum -y install python-pip
RUN yum -y install java-1.8.0-openjdk
RUN yum -y install perl

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /app
RUN wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip && unzip bcl2fastq*.zip && rpm -i *.rpm
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && unzip fastqc*.zip

ADD . /app

RUN chmod 755 /app/FastQC/fastqc
RUN ln -s /app/FastQC/fastqc /usr/local/bin/fastqc

RUN pip install --upgrade pip
RUN pip install numpy jinja2 MySQL-python

ENV PYTHONPATH=/app

CMD ["python", "odybcl2fastq/process_runs.py"]
