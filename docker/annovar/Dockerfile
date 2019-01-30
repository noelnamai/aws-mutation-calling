#specifying the base image
FROM ubuntu:latest

#name and email of the person who maintains the file
MAINTAINER Noel Namai namai@bcm.edu

#Updating ubuntu and installing other necessary software
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

#install ANNOVAR
COPY annovar.tar.gz .
RUN tar -xzf annovar.tar.gz \
	&& rm annovar.tar.gz \
	&& mv annovar /opt/

#chmod of /opt/
RUN chmod -R 777 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
