#######################################################################
#     Basic image
#######################################################################
FROM cgap/cgap-ubuntu2004-py-38:0.0.1
MAINTAINER Phil Grayson (phil.d.grayson@gmail.com), Michele Berselli (berselli.michele@gmail.com)

#######################################################################
#     General updates & installing necessary Linux components
#     VEP:
#       - libbz2-dev, liblzma-dev (HTSlib)
#       - libdbi-perl, libmodule-build-perl,
#         libarchive-zip-perl, libtest-exception-perl (Perl)
#       - libmysqlclient-dev, libpng-dev,
#         libexpat1-dev, libssl-dev (other)
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix \
    libbz2-dev \
    liblzma-dev \
    libdbi-perl \
    libmodule-build-perl \
    libarchive-zip-perl \
    libtest-exception-perl \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    libexpat1-dev

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

#######################################################################
#     Software
#######################################################################
## VEP
RUN git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    git checkout release/101 && \
    perl INSTALL.pl --AUTO a --NO_UPDATE && \
    cd ..

## perl module for VEP
RUN cpan App::cpanminus && \
    wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz && \
    tar xzf v335_base.tar.gz && \
    export KENT_SRC=/usr/local/bin/kent-335_base/src && \
    export MACHTYPE=$(uname -m) && \
    export CFLAGS="-fPIC" && \
    MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'` && \
    export MYSQLINC && \
    MYSQLLIBS=`mysql_config --libs` && \
    export MYSQLLIBS && \
    cd $KENT_SRC/lib && \
    echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk && \
    make clean && make && \
    cd ../jkOwnLib && \
    make clean && make && \
    cd /usr/local/bin/ && \
    mkdir -p /usr/local/bin/cpanm_bw && \
    export PERL5LIB=$PERL5LIB:/usr/local/bin/cpanm_bw/lib/perl5 && \
    cpanm -l /usr/local/bin/cpanm_bw Bio::DB::BigFile

## Test
RUN export PERL5LIB=$PERL5LIB:/usr/local/bin/ensembl-vep && \
    cd /usr/local/bin/ensembl-vep && \
    perl -Imodules t/AnnotationSource_File_BigWig.t && \
    cd /usr/local/bin

#######################################################################
#     Scripts
#######################################################################
## vep-annot
COPY scripts/vep-annot.sh .
RUN chmod +x vep-annot.sh

#######################################################################
#     Setting env variables
#######################################################################
ENV PATH=/usr/local/bin/ensembl-vep/:$PATH
ENV PERL5LIB=$PERL5LIB:/usr/local/bin/cpanm_bw/lib/perl5:/usr/local/bin/ensembl-vep

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]