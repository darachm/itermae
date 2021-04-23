FROM python:3.9-buster AS itermae

## Setup, i8n and base packages
#RUN apt -y update
#RUN apt-get -y install apt-utils software-properties-common language-pack-en wget make
#RUN update-locale

LABEL version=0.6.0.1
ENV DEBIAN_FRONTEND="noninteractive"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV LANGUAGE=C.UTF-8
RUN git clone https://gitlab.com/darachm/itermae.git 
WORKDIR /itermae
RUN make install
WORKDIR /
RUN rm -rf itermae

FROM itermae AS itermae-plus
WORKDIR /
RUN apt-get -y install wget 
RUN wget https://ftp.gnu.org/gnu/parallel/parallel-20201222.tar.bz2
RUN tar -xvjf parallel-20201222.tar.bz2
WORKDIR /parallel-20201222 
RUN ./configure 
RUN make 
RUN make install 
RUN make clean
RUN mkdir /root/.parallel
RUN touch /root/.parallel/will-cite
WORKDIR /
RUN rm -rf /parallel-20201222
RUN rm -rf /parallel-20201222.tar.bz2
RUN apt-get -y install gzip xz-utils bzip2 unzip curl libcurl4-openssl-dev mawk perl

