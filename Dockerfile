FROM ubuntu:22.04
LABEL Description="rdkit build environment"

ENV HOME /root

RUN apt-get update && apt-get -y install cmake && apt-get -y install build-essential
RUN apt-get -y install libasio-dev
RUN apt-get -y install libgflags-dev
RUN apt-get install -y libfreetype-dev libfreetype6 libfreetype6-dev libfontconfig1-dev
RUN apt-get -y install curl
RUN apt-get -y install git

RUN curl -O https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
RUN chmod +x Anaconda3-2024.02-1-Linux-x86_64.sh
RUN ./Anaconda3-2024.02-1-Linux-x86_64.sh -b

ENV PATH /root/anaconda3/bin:$PATH

RUN conda update -y conda
RUN conda update -y --all

RUN conda create -y -n my-rdkit-env
RUN conda init

#RUN conda activate my-rdkit-env
#RUN conda install -y numpy 
#RUN conda install -y cmake cairo pillow eigen pkg-config
#RUN conda install -y boost-cpp boost py-boost
#RUN conda install -y -c conda-forge libstdcxx-ng ????
#RUN conda install -c conda-forge libstdcxx-ng=12 - works

ENV RDBASE /src/cpp/third_party/rdkit

WORKDIR /src/cpp/third_party/rdkit

RUN conda install -y numpy 
RUN conda install -y cmake cairo pillow eigen pkg-config
RUN conda install -y boost-cpp boost py-boost
RUN conda install -y -c conda-forge libstdcxx-ng 
#COPY conda_shell.sh .
#RUN ./conda_shell.sh

SHELL ["/bin/bash", "-c"]