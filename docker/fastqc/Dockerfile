#specifying the base image
FROM java:8

#name and email of the person who maintains the file
MAINTAINER Noel Namai noelnamai@yahoo.com

# Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
	&& apt-get clean

#install Fastqc-0.11.3
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip \
	&& unzip fastqc_v0.11.3.zip \
	&& mv FastQC/* /usr/bin/ \
	&& rm -rf FastQC fastqc_v0.11.3.zip

#chmod of /usr/bin/ and change working directory
RUN chmod -R 777 /usr/bin/
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
