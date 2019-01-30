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

#install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
	&& unzip Trimmomatic-0.38.zip \
	&& mv Trimmomatic-0.38/* /opt/ \
	&& rm Trimmomatic-0.38.zip

#chmod of /opt/ and change working directory
RUN chmod -R 777 /opt/
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
