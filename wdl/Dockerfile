FROM ubuntu:20.04

RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev libboost-all-dev  && \
	apt-get clean



WORKDIR /opt


### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
  bunzip2 htslib-1.17.tar.bz2 && \
  tar xvf htslib-1.17.tar && \
	cd htslib-1.17/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make -j 4 && \
	make install 
	
	 

#	cp lib/* /usr/lib/ && \

#3cp -r include/* /usr/include/

ADD https://api.github.com/repos/chaissonlab/vamos/git/refs/heads/ version.json

RUN wget https://github.com/ChaissonLab/countme/archive/refs/tags/countme-v0.1.tar.gz
RUN tar zxvf  countme-v0.1.tar.gz
RUN ls 
RUN cd countme-countme-v0.1/ && make CONDA_PREFIX=/usr/local
RUN cd countme-countme-v0.1/ && cp countme /usr/bin

ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"
