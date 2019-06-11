#
# Docker container for dev odybcl2fastq
#
FROM conda/miniconda3-centos7

RUN yum -y update
RUN yum install -y gcc
RUN yum install -y sssd

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /app

RUN conda install -y jinja2 numpy
RUN conda install -y -c bioconda mysqlclient
RUN conda install -y -c bioconda -c conda-forge snakemake

ENV PYTHONPATH=/app
ADD . /app
RUN chown -R 559147:403265 /app
WORKDIR /app/odybcl2fastq

CMD ["python", "/app/odybcl2fastq/process_snakemake_runs.py"]
