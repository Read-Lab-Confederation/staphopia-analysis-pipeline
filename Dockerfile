FROM broadinstitute/gatk3:3.8-0 as gatk
FROM rpetit3/nextconda-base
MAINTAINER robbie.petit@gmail.com

# Install dependencies
# Install via Bioconda
# Program versions last update 11/2017
# Update Nextflow and Conda
RUN apt-get -qq update \
    && apt-get -qq -y --no-install-recommends install \
        gcc \
        g++ \
        less \
        libicu-dev \
        libxml2-dev \
        python3-tk \
        wget \
        zlib1g-dev \
    && conda upgrade conda \
    && conda install -y nextflow=0.28.2 \
    && conda install -y bbmap=37.66 \
    && conda install -y bedtools=2.26.0gx \
    && conda install -y bwa=0.7.17 \
    && conda install -y jellyfish=2.2.6 \
    && conda install -y picard=2.14.1 \
    && conda install -y samtools=1.6 \
    && conda install -y spades=3.11.1 \
    && conda install -y mentalist=0.1.3 \
    # ARIBA
    && conda install -y bowtie2=2.3.3.1 \
    && conda install -y cd-hit=4.6.8 \
    && conda install -y mummer=3.23 \
    && conda install -y 'icu=56.*' \
    # PROKKA
    && conda install -y perl-bioperl=1.6.924 \
    && conda install -y perl-xml-simple \
    && conda install -y tbl2asn=25.6 \
    # Other
    && conda install -y fastq-scan=0.3 \
    && conda install -y assembly-scan=0.2 \
    && conda install -y vcf-annotator=0.5 \
    && conda clean --all --yes \
    && pip install --upgrade pip setuptools\
    && pip install ariba \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log /tmp/* /var/tmp/* \
    && mkdir -p /data /tmp/install

# Final Programs
# Aspera Connect
RUN cd /tmp/install \
    && wget --quiet http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz \
    && tar -xzf aspera-connect-3.7.4.147727-linux-64.tar.gz \
    && sed -i 's=INSTALL_DIR\=~/.aspera/connect=INSTALL_DIR\=/opt/aspera=' aspera-connect-3.7.4.147727-linux-64.sh \
    && mkdir -p /opt/aspera \
    && bash aspera-connect-3.7.4.147727-linux-64.sh \
# BLAST 2.7.1
    && cd /tmp/install \
    && curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz -o ncbi-blast-2.7.1+.tar.gz \
    && tar -xzf ncbi-blast-2.7.1+.tar.gz \
    && cd /tmp/install/ncbi-blast-2.7.1+/bin \
    && mv blastn blastp makeblastdb tblastn /usr/local/bin \
# ENA Downloader (ena-dl)
    && cd /tmp/install \
    && curl -sSL https://github.com/rpetit3/ena-dl/archive/v0.1.tar.gz -o ena-dl-0.1.tar.gz \
    && tar -xzf ena-dl-0.1.tar.gz \
    && pip install -r ena-dl-0.1/requirements.txt \
    && chmod 755 ena-dl-0.1/ena-* \
    && mv ena-dl-0.1/ena-* /usr/local/bin/ \
# Illumina Cleanup
    && cd /tmp/install \
    && curl -sSL https://github.com/rpetit3/illumina-cleanup/archive/v0.3.tar.gz -o illumina-cleanup-0.3.tar.gz \
    && tar -xzf illumina-cleanup-0.3.tar.gz \
    && pip install -r illumina-cleanup-0.3/requirements.txt \
    && chmod 755 illumina-cleanup-0.3/src/* \
    && mv illumina-cleanup-0.3/src/*.py /usr/local/bin/ \
    && g++ -Wall -O3 -o /usr/local/bin/fastq-interleave illumina-cleanup-0.3/src/fastq-interleave.cpp \
    && mkdir -p /opt/staphopia/data/fastq \
    && mv illumina-cleanup-0.3/data/*.fasta /opt/staphopia/data/fastq \
# PROKKA
    && cd /tmp/install \
    && curl -sSL https://github.com/rpetit3/prokka/archive/v1.12-staphopia.tar.gz -o prokka-1.12-staphopia.tar.gz \
    && tar -xzf prokka-1.12-staphopia.tar.gz \
    && mv prokka-1.12-staphopia/ /opt/prokka \
    && export PATH=/opt/prokka/bin:$PATH

ENV ASCP /opt/aspera/bin/ascp
ENV ASCP_KEY /opt/aspera/etc/asperaweb_id_dsa.openssh
ENV PATH $PATH:/opt/prokka/bin

# Final touches
# GATK
COPY --from=gatk /usr/GenomeAnalysisTK.jar /usr/local/bin/GenomeAnalysisTK.jar
COPY data /tmp/data
RUN mkdir -p /opt/staphopia/data \
    && cd /tmp/data/ \
    && ls *.tar.gz | xargs -I {} tar xzf {} \
    && rm *.tar.gz \
    && mv ./* /opt/staphopia/data/ \
    && prokka --setupdb

COPY scripts /tmp/scripts
RUN chmod 755 /tmp/scripts/* \
    && mv /tmp/scripts/* /usr/local/bin \
    && rm -rf /tmp/*

# tbl2asn Updates, last updated 20190306
run cd /tmp/ \
    && wget --quiet -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip tbl2asn.gz \
    && chmod 755 tbl2asn \
    && mv tbl2asn /usr/local/bin/tbl2asn \
    && ln -f -s /usr/local/bin/tbl2asn /opt/prokka/binaries/linux/tbl2asn

WORKDIR /data
