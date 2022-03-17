FROM 4dndcic/ubuntu16.04-miniconda-python36:v1
MAINTAINER Phil Grayson (phil.d.grayson@gmail.com), Michele Berselli (berselli.michele@gmail.com)

# general updates & installing necessary Linux components
RUN apt-get update && apt-get install -y \
    nano \
    htop \
    bzip2 \
    gcc \
    g++ \
    git \
    less \
    libncurses-dev \
    make \
    time \
    zip \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    liblz4-tool \
    libbz2-dev \
    liblzma-dev \
    locales \
    parallel \
    autoconf \
    autogen \
    autotools-dev \
    automake \
    pkg-config \
    tabix \
    yasm \
    libdbi-perl \
    libmodule-build-perl \
    libarchive-zip-perl \
    libmysqlclient-dev \
    libpng-dev \
    libexpat1-dev

WORKDIR /usr/local/bin

# vcftools for vcf-integrity-check script
RUN git clone https://github.com/vcftools/vcftools.git && \
    mv vcftools vcftools_tmp && \
    cd vcftools_tmp && \
    git checkout 954e607 && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cd ..

# multiple_vcf-integrity-check script
COPY multiple_vcf-integrity-check.sh .
RUN chmod +x multiple_vcf-integrity-check.sh

# python and pysam (required for VEP)
RUN conda install python=3.6.8 && \
    conda clean -a -y -f

# granite
RUN pip install granite-suite==0.1.13

# snv_somatic scripts: somatic_sentieon_vcf_splitter.py SV_bnd_converter.py
COPY scripts/somatic_sentieon_vcf_splitter.py .
RUN chmod +x somatic_sentieon_vcf_splitter.py
COPY scripts/SV_bnd_converter.py .
RUN chmod +x SV_bnd_converter.py

# supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
