#######################################################################
#     Basic image
#######################################################################
FROM cgap/cgap-ubuntu2004-py-38:0.0.1
MAINTAINER Michele Berselli (berselli.michele@gmail.com), Dominika Maziec (dominika_maziec@hms.harvard.edu)

#######################################################################
#     General updates & installing necessary Linux components
#       - autoconf, pkg-config
#.      - bcftools
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix \
    autoconf \
    pkg-config \
    bcftools

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

#######################################################################
#     Software
#######################################################################
## granite pandas
RUN pip install granite-suite==0.2.0 pandas==1.5.2

#######################################################################
#     Scripts
#######################################################################

## snv_somatic scripts: somatic_annot.py
COPY scripts/somatic_annot.py .
RUN chmod +x somatic_annot.py

#######################################################################
#     Setting env variables
#######################################################################
## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
