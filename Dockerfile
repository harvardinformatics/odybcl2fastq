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
RUN wget -O bcl2fastq2-v2-20-0-linux-x86-64.zip "https://files.softwaredownloads.illumina.com/e8ed3335-5201-48ff-a2bc-db4bfb792c85/bcl2fastq2-v2-20-0-linux-x86-64.zip?Expires=1573826916&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9maWxlcy5zb2Z0d2FyZWRvd25sb2Fkcy5pbGx1bWluYS5jb20vZThlZDMzMzUtNTIwMS00OGZmLWEyYmMtZGI0YmZiNzkyYzg1L2JjbDJmYXN0cTItdjItMjAtMC1saW51eC14ODYtNjQuemlwIiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNTczODI2OTE2fX19XX0_&Signature=A19lLm-kA2NYGYmUkc1l-b39UrnIzL525SVA~6frAris4gi-S~mGt8C4HhkW9komKivqYlD3VrqdTKmkRzOS9CRD0F5xrFLlJ3~D37U~KlHQy0aCXsFW64~zl-dIDZ7FkojK1gxNGKfWlgs8xewUl75oVie3bPrIavqYzWVRq9XNWk-z0WUosCLa1ih7SIXX7I0uogvR7AAFq30FzI5~EGidOBh88o5Mftt3AfI8vfPeFw1P6UW4~Fa4fr0kBRcvKiTx0FmtOYV5hzfm-q-jFSCuKU8BdrXLdxfmmMYWlx7tBOjrUgqCt9~Krdi59~w5f8GmLMCVPW4zBIbpEz8UIQ__&Key-Pair-Id=APKAJO3UYWPXK4A26FPQ" && unzip bcl2fastq*.zip && rpm -i *.rpm
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && unzip fastqc*.zip

ADD . /app

RUN chmod 755 /app/FastQC/fastqc
RUN ln -s /app/FastQC/fastqc /usr/local/bin/fastqc

RUN yum -y install vim
RUN pip install --upgrade pip
RUN pip install numpy jinja2 MySQL-python

ENV PYTHONPATH=/app

CMD ["python", "odybcl2fastq/process_runs.py"]
