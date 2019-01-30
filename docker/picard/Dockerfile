#specifying the base image
FROM java:8

#name and email of the person who maintains the file
MAINTAINER Noel Namai noelnamai@yahoo.com

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

#install picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.18.25/picard.jar \
	&& mv picard.jar /opt/

#chmod of /opt/
RUN chmod -R 777 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
