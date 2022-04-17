ARG VERSION=20.04
FROM ubuntu:$VERSION

MAINTAINER eab <artscan@list.ru>

RUN echo 'APT::Get::Assume-Yes "true";' >> /etc/apt/apt.conf \
    && apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install \
    linux-headers-$(uname -r) \
    zip \
    git \
    mlocate \
    sudo \
    cmake \
    wget \
    maxima \
    libblas-dev \
    liblapack-dev \
    libblas3 \
    liblapack3 \
    libopenblas-dev \
    libopenblas-base \
    openmpi-bin \
    mpich \
    libopenmpi-dev \
    libmetis-dev \
    gfortran \
    python-dev \
    python-numpy \
    libreadline5 \
    libffcall1-dev \
    libavcall1 \
    libatlas-base-dev \
    libsigsegv-dev \
    && wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python get-pip.py \
    # su-exec
    && git clone https://github.com/ncopa/su-exec.git /tmp/su-exec \
    && cd /tmp/su-exec \
    && make \
    && chmod 770 su-exec \
    && mv ./su-exec /usr/local/sbin/ \
    && apt-get autoremove \
    && rm -rf /tmp/* /var/lib/apt/lists/* /root/.cache/*

RUN pip2 install scipy

COPY asEnvUser /usr/local/sbin/
# Only for sudoers
RUN chown root /usr/local/sbin/asEnvUser \
    && chmod 700  /usr/local/sbin/asEnvUser

RUN updatedb
RUN useradd -u 1000 eab
RUN groupmod -g 1000 eab

COPY . /home/eab/git/difwave/abiem/.
RUN chown -R eab:eab /home/eab/git/difwave/abiem/

ENV UNAME="eab" \
    GNAME="eab" \
    UHOME="/home/eab" \
    UID="1000" \
    GID="1000" \
    SHELL="/bin/bash" \
    WORKSPACE="/home/eab/git/difwave/abiem"

WORKDIR "${WORKSPACE}"

ENTRYPOINT ["asEnvUser"]


