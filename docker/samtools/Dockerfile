#specifying the base image
FROM ubuntu:16.04

#name and email of the person who maintains the file
MAINTAINER Noel Namai namai@bcm.edu

# Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		git \
		g++ \
		autoconf \
		automake \
		gcc \
		perl \
		libbz2-dev \
		liblzma-dev \
		libcurl4-gnutls-dev \
		libssl-dev \
		libncurses5-dev \
		make \
		ncurses-dev \
		pkg-config \
		zlib1g-dev \
		cmake \
	&& apt-get clean

#install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
	&& tar -xvjf htslib-1.9.tar.bz2 \
	&& rm htslib-1.9.tar.bz2 \
	&& cd htslib-1.9 \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#install Samtools
RUN wget https://github.com/samtools/samtools/archive/1.9.tar.gz \
	&& tar -xvzf 1.9.tar.gz \
	&& rm 1.9.tar.gz \
	&& cd samtools-1.9 \
	&& autoheader \
	&& autoconf -Wno-syntax \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#install VCFtools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
	&& tar -xvzf vcftools-0.1.16.tar.gz \
	&& rm vcftools-0.1.16.tar.gz \
	&& cd vcftools-0.1.16 \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#chmod of /opt/
RUN chmod -R 777 /opt/

#change working directory to /opt/
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
