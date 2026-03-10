FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Allow OpenMPI to run inside Docker
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Install basic tools and MPI/NetCDF stack
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        environment-modules \
        cmake \
        pkg-config \
        gcc g++ gfortran \
        python3 python3-pip python3-venv \
        make git wget \
        openmpi-bin openmpi-common libopenmpi-dev \
        libnetcdf-dev netcdf-bin \
        libnetcdff-dev libnetcdff7 \
        libpnetcdf-dev pnetcdf-bin \
        vim ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# aliases 
RUN mkdir -p /usr/local/bin && \
    ln -s $(which bash) /usr/local/bin/sbatch && \
    ln -s $(which python3) /usr/bin/python

# Trick miniweather scripts into thinking Slurm is available
ENV SLURM_JOB_ID=1

# Prevent permission issues in host by creating user in container
ARG UID=1000
ARG GID=1000

RUN groupadd -g $GID user && \
    useradd -m -u $UID -g $GID -s /bin/bash user

USER user
WORKDIR /workspace

CMD ["/bin/bash"]
