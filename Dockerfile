#
# Docker container for dev odybcl2fastq
#
FROM centos:7.7.1908

RUN yum -y install \
  MySQL-python \
  java-1.8.0-openjdk \
  numpy \
  perl \
  python \
  python-jinja2 \
  unzip \
  vim \
  wget

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /app
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && unzip fastqc*.zip && rm -f fastqc*.zip
RUN chmod 755 /app/FastQC/fastqc
RUN ln -s /app/FastQC/fastqc /usr/local/bin/fastqc

ADD . /app

# assumes odybcl2fastq/config.json in build context
RUN [ -s odybcl2fastq/config.json ]
# assumes bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm is in build context
RUN rpm -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm -f bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
RUN ln -s /tmp/centrifuge.log /app/centrifuge.log 

ENV PYTHONPATH=/app

CMD ["python", "odybcl2fastq/process_runs.py"]
