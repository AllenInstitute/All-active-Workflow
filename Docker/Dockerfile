# Docker Image for Optimization

FROM continuumio/miniconda3
MAINTAINER Ani Nandi <anin@alleninstitute.org>


RUN apt-get update 
RUN apt-get install -y wget libx11-6 python-dev git build-essential         \
    libncurses-dev autoconf automake gcc g++ make gfortran nano

# Install MPICH

RUN wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz              

RUN                                                                         \
    tar xzf mpich-3.2.tar.gz                                            &&  \
    cd mpich-3.2                                                        &&  \
    ./configure                                                         &&  \
    make -j 4                                                           &&  \
    make install                                                        
RUN rm -rf mpich-3.2 && rm mpich-3.2.tar.gz


# # ENV _USER=docker
# # ENV _USER_ID=1000
ENV HOME=/home

# RUN useradd -m $_USER && echo $_USER:$_USER | chpasswd && adduser $_USER sudo
# USER $_USER
#RUN chown -R $_USER_ID $HOME

ENV BUILD_DIR=$HOME/build
RUN mkdir -p $BUILD_DIR 

WORKDIR $BUILD_DIR


# Python packages

RUN conda create -n aibs_opt python=3.7 pytables
RUN echo "source activate aibs_opt" > ~/.bashrc

ENV CONDA_ENV_PATH=/opt/conda/envs/aibs_opt
ENV PATH=$CONDA_ENV_PATH/bin:$PATH


# Neuron 7.6

ENV NRN_VER=7.6
ENV NRN=nrn-$NRN_VER
ENV NRN_INSTALL_DIR=$BUILD_DIR/$NRN

WORKDIR $BUILD_DIR

RUN wget http://www.neuron.yale.edu/ftp/neuron/versions/v$NRN_VER/$NRN.tar.gz; \
	tar xzf $NRN.tar.gz && rm $NRN.tar.gz

WORKDIR $NRN_INSTALL_DIR

RUN ./configure --with-nrnpython=$CONDA_ENV_PATH/bin/python --with-paranrn --without-iv --prefix=$(pwd); \
    make -j 4; make install

RUN cd src/nrnpython; python setup.py install
RUN cd $CONDA_ENV_PATH/bin; ln -s $NRN_INSTALL_DIR/x86_64/bin/nrnivmodl;\
    ln -s $NRN_INSTALL_DIR/x86_64/bin/nrniv

COPY . /app
WORKDIR /app
RUN pip install --upgrade pip                                               

RUN pip install -r requirements.txt

# Install mpi4py

RUN wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.1.tar.gz       

RUN                                                                         \
    tar zxf mpi4py-3.0.1.tar.gz                                         &&  \
    cd mpi4py-3.0.1                                                     &&  \
    python setup.py build                                               &&  \
    python setup.py install                                             

RUN rm -rf mpi4py-3.0.1

#WORKDIR $HOME

