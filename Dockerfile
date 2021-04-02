FROM ubuntu:focal

RUN export LC_ALL=C.UTF-8
RUN export LANG=C.UTF-8
RUN export DEBIAN_FRONTEND="noninteractive"

RUN apt -y update
RUN apt -y install apt-utils software-properties-common language-pack-en \
        wget make
RUN update-locale

RUN wget https://ftp.gnu.org/gnu/parallel/parallel-20201222.tar.bz2
RUN tar -xvjf parallel-20201222.tar.bz2
RUN cd parallel-20201222
RUN ./configure && make && make install
RUN cd /
RUN echo "will cite\n" | parallel --citation || echo "why is this a return value of 1"

RUN apt -y install python3 python3-pip
RUN apt -y install gzip xz-utils bzip2 zip curl libcurl4-openssl-dev gawk perl

RUN python3 -m pip install --upgrade pip

LABEL version=0.5.1
RUN python3 -m pip install itermae==0.5.1

CMD itermae
