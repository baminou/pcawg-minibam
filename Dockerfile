FROM ubuntu:16.04
LABEL MAINTAINER solomon.shorser@oicr.on.ca
RUN mkdir -p /opt/gitroot

RUN apt-get update
#RUN apt list upgradable
#RUN apt-get upgrade -y

RUN apt-get install -y apt-utils git tabix bash=4.3-14ubuntu1.2

# Install bcftools - this is used for normalizing VCFs

RUN cd /opt/gitroot && git clone https://github.com/samtools/htslib.git

RUN cd /opt/gitroot && git clone https://github.com/samtools/bcftools.git

RUN cd /opt/gitroot && git clone https://github.com/vcftools/vcftools.git

RUN apt-get install -y gcc make zlib1g-dev libstring-random-perl

RUN cd /opt/gitroot/htslib \
	&& git checkout 1.3.1 \
	&& make && make install

RUN cd /opt/gitroot/bcftools \
	&& git checkout 1.3.1 \
	&& make && make install

# Instal VCFTools - this is used for vcf-sort in the vcf_merge_by_type.pl script
RUN apt-get install pkg-config dh-autoreconf -y

RUN cd /opt/gitroot/vcftools \
	&& git checkout v0.1.14 \
	&& ./autogen.sh && ./configure && make && make install

# Copy in the scripts.
RUN mkdir -p /opt/oxog_scripts
COPY ./scripts/* /opt/oxog_scripts/
RUN chmod a+x /opt/oxog_scripts/*.*
