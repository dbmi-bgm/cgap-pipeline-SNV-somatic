FROM cgap/cgap-ubuntu2004-py-38:0.0.1
MAINTAINER Phil Grayson (phil.d.grayson@gmail.com), Michele Berselli (berselli.michele@gmail.com)

RUN apt-get update -y && apt-get install -y \
    autoconf \
    pkg-config \
    locales \
    tabix

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
