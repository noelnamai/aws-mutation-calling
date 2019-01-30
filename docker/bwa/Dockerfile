#specifying the base image
FROM ubuntu:16.04

#name and email of the person who maintains the file
MAINTAINER Noel Namai noelnamai@yahoo.com

# Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get install -y \
		bwa \
	&& apt-get clean


#chmod of /opt/ and change working directory
RUN chmod -R 777 /opt/
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
