# Development image for Process, used in challenging environments to give
# identical results to CI or to quickly provide a working environment to a
# casual user. Differs from CI image by not using virtual environments for
# differt Python versions and doesn't require the CI script to perform extra
# installations
FROM ubuntu:20.04
LABEL description="Process development image"

# Required to get around interactive timezone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV FC=gfortran

# Install apt packages
# poppler-utils is required by pdf2image
RUN apt-get update && apt-get install -y \
    gdb \
    gfortran \
    git \
    graphviz \
    libblas-dev \
    liblapack-dev \
    liblapack3 \
    poppler-utils \
    python3-pip \
    python-is-python3 \
    libhdf5-serial-dev \
    netcdf-bin \
    libnetcdf-dev

# Install Python packages into a virtual environment
# Use virtualenv to prevent "running pip as root" warnings
RUN python -m pip install --upgrade pip && \
    python -m pip install \
    cmake \
    Markdown \
    matplotlib \
    mkdocs \
    mkdocs-material \
    netCDF4 \
    numpy \
    pathlib \
    pillow \
    pymdown-extensions \
    pytest \
    scipy \
    toposort \
    wheel
