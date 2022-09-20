# Set the base image to Ubuntu
FROM ubuntu:16.04

RUN apt-get update && \
    apt-get install -y redis-server && \
    apt-get clean

EXPOSE 6379

CMD ["redis-server", "--protected-mode no"]

# File Author / Maintainer
MAINTAINER jvn2 <noort@physics.leidenuniv.nl>

ARG GUPPY_VERSION=5.0.7
# ARG MEGALODON_VERSION=2.3.5
ARG MEGALODON_VERSION=2.4.1

#ARG CONDA_VERSION=py38_4.9.2
ARG CONDA_MD5=122c8c9beb51e124ab32a0fa6426c656

ARG CONDA_VERSION=py37_4.11.0

ARG SAMTOOLS_VERSION=1.11

ARG TOMBO_VERSION=1.5.1

WORKDIR /home

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update -y && \
    apt-get install --yes \
                        wget \
                        apt-transport-https \
                        bzip2 \
                        ca-certificates \
                        libglib2.0-0 \
                        libsm6 \
                        libxext6 \
                        libxrender1 \
                        mercurial \
                        subversion \
                        libcurl4-openssl-dev \
                        libssl-dev \
                        libhdf5-cpp-11 \
                        libzmq5 \
                        libboost-atomic1.58.0 \
                        libboost-chrono1.58.0 \
                        libboost-date-time1.58.0 \
                        libboost-filesystem1.58.0 \
                        libboost-program-options1.58.0 \
                        libboost-regex1.58.0 \
                        libboost-system1.58.0 \
                        libboost-log1.58.0 \
                        libboost-iostreams1.58.0 \
                        build-essential \
                        git \
                        libgl1-mesa-glx


#RUN /bin/sh -c apt-get update && apt-get install -y     build-essential     zlib1g     zlib1g-dev     git     perl     python-pip     make     vim     wget     gcc     g++     perl     zip     gzip     gawk # buildkit
#RUN /bin/sh -c apt-get update && apt-get install -y bzip2 wget libgl1-mesa-glx git
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh
RUN bash Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -b -p miniconda3 # buildkit
ENV PATH=/usr/local/bin/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV HDF5_PLUGIN_PATH=/usr/local/bin/miniconda3/hdf5/lib/plugin

RUN /bin/sh -c conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels anaconda # buildkit
RUN /bin/sh -c conda install -c bioconda ont-tombo # buildkit
RUN /bin/sh -c conda install -c bioconda ont_vbz_hdf_plugin # buildkit
RUN /bin/sh -c conda install -c bioconda meme # buildkit
RUN /bin/sh -c pip install "numpy<1.20" # buildkit
RUN /bin/sh -c pip install ont-fast5-api # buildkit
RUN /bin/sh -c ln -s /usr/local/bin/miniconda3/bin/* /usr/local/bin # buildkit

RUN cd /tmp &&\
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_cpu_${GUPPY_VERSION}-1~xenial_amd64.deb && \
    dpkg -i /tmp/ont_guppy_cpu_${GUPPY_VERSION}-1~xenial_amd64.deb && \
    rm *.deb && \

    #Newest version of pysam does not work with bzip2
    pip3 install pysam==0.18.0 && \
    pip3 install numpy cython ont_pyguppy_client_lib==${GUPPY_VERSION} && \
    pip3 install megalodon==${MEGALODON_VERSION} && \
    pip3 install ont-tombo==${TOMBO_VERSION} && \

    #install samtools so megalodon works with --sort-mappings
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --without-curses --disable-bz2 --disable-lzma && \
    make  && \
    make install  && \
    cd ..  && \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2 && \

    git clone https://github.com/nanoporetech/rerio /home/rerio && \
    /home/rerio/download_model.py /home/rerio/basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001 && \
    /home/rerio/download_model.py /home/rerio/basecall_models/res_dna_r941_min_modbases_5mC_5hmC_CpG_v001 && \
    /home/rerio/download_model.py /home/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 && \

    mkdir -p /home/ont-guppy/bin && \
    ln -s /usr/bin/guppy_basecall_server /home/ont-guppy/bin/guppy_basecall_server && \

    mkdir /home/nanopore && \
    git clone https://github.com/JvN2/Nanopore.git /home/nanopore && \

    mkdir /home/data && \

    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*



CMD ["echo", "Image created"]