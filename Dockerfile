FROM ubuntu:focal

# Setup, i8n and base packages
RUN export LC_ALL=C.UTF-8
RUN export LANG=C.UTF-8
RUN export DEBIAN_FRONTEND="noninteractive"
RUN apt -y update
RUN apt -y install apt-utils software-properties-common language-pack-en \
        wget make
RUN update-locale

# Install GNU parallel from source because pkg manager is odd
RUN wget https://ftp.gnu.org/gnu/parallel/parallel-20201222.tar.bz2
RUN tar -xvjf parallel-20201222.tar.bz2
RUN cd parallel-20201222 && ./configure && make && make install
RUN echo "will cite\n" | parallel --citation || echo "why is this a return value of 1"

# Installing python, pip, some utilities
RUN apt -y install python3 python3-pip
RUN apt -y install gzip xz-utils bzip2 zip curl libcurl4-openssl-dev gawk perl
RUN python3 -m pip install --upgrade pip
RUN apt -y install git

LABEL version=0.5.1
RUN git clone https://gitlab.com/darachm/itermae.git 
RUN cd itermae && make install #python3 -m pip install itermae==0.5.1

CMD itermae
