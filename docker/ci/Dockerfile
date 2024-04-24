# Docker image for CI system
FROM ubuntu:20.04
LABEL description="Process Ubuntu 20.04 image for running CI jobs"


# Required to get around interactive timezone prompt
ENV DEBIAN_FRONTEND=noninteractive


COPY requirements.txt /tmp/requirements.txt
COPY requirements_dev.txt /tmp/requirements_dev.txt


# Install apt packages
# gfortran, debugging tools, git, maths and Latex libraries
# poppler-utils is required by pdf2image
# This syntax is recommended by Docker to bust caching
RUN apt-get update && apt-get install -y \
    gdb \
    gfortran \
    git \
    graphviz \
    libblas-dev \
    liblapack-dev \
    liblapack3 \
    pandoc \
    poppler-utils \
    python3-pip \
    python-is-python3 \
    texlive-base \
    texlive-font-utils \
    texlive-latex-extra \
    lcov \
    cmake \
    bc


# Install apt packages for pyenv.
RUN apt-get install -y \
    build-essential \
    curl \
    libbz2-dev \
    libffi-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    llvm \
    python-openssl \
    tk-dev \
    wget \
    xz-utils \
    zlib1g-dev


# Install pyenv.
RUN curl https://pyenv.run | bash
ENV PYENV_ROOT /root/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PYENV_ROOT/versions/3.10.14/bin:$PATH

RUN pyenv install 3.10.14 && \
    pyenv global 3.10.14


RUN python -m pip install --upgrade pip && \
    python -m pip install -r /tmp/requirements.txt && \
    python -m pip install -r /tmp/requirements_dev.txt && \
    python -m pip install git+https://github.com/jonmaddock/ford

RUN git clone https://github.com/jmorris-uk/anybadge.git /root/anybadge/ && pip install anybadge
