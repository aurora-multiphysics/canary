# syntax=docker/dockerfile:1
#
# canary-moose-dev
#
# Everything needed to build MOOSE with MFEM support: PETSc, libMesh, WASP,
# Conduit, MFEM and the MOOSE python environment, each at the version pinned by
# MOOSE_REVISION. The steps follow MOOSE's own apptainer/{mpi,petsc,libmesh,
# moose-dev}.def, but start from IDAHO's published compiler and MPI image
# rather than the intermediate images that MOOSE only publishes to INL's
# registry. Only the 'opt' method is built.
#
# Required build arguments:
#
#   MOOSE_REVISION: MOOSE commit whose submodules and build scripts are used
#

ARG MOOSE_MPI_IMAGE=ghcr.io/idaholab/moose-containers/moose-mpi-ubuntu24-gcc:24.04-gcc14.2.0-mpich5.0.2-openmpi5.0.11-20261008

# MOOSE source and the dependency submodules. The build steps below mount this
# read-write; anything they write to it is discarded with the step, so none of
# the source or build trees end up in the image.
# hadolint ignore=DL3006
FROM ${MOOSE_MPI_IMAGE} AS source
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG MOOSE_REVISION
# MOOSE marks the mfem and conduit submodules 'update = none', which a plain
# 'git submodule update' silently skips, so they need --checkout. Each
# submodule is then checked for content so that one that is skipped fails here
# rather than in its build step.
RUN test -n "$MOOSE_REVISION" && \
    git clone --filter=blob:none --no-checkout https://github.com/idaholab/moose.git /moose && \
    git -C /moose checkout --quiet "$MOOSE_REVISION" && \
    git -C /moose submodule update --init --depth 1 petsc && \
    git -C /moose submodule update --init --recursive --depth 1 \
        libmesh framework/contrib/wasp && \
    git -C /moose submodule update --init --checkout --depth 1 framework/contrib/mfem && \
    git -C /moose submodule update --init --checkout --recursive --depth 1 \
        framework/contrib/conduit && \
    for submodule in petsc libmesh framework/contrib/wasp framework/contrib/mfem \
        framework/contrib/conduit; do \
        test -n "$(ls -A "/moose/${submodule}")" || \
            { echo "Submodule ${submodule} was not checked out"; exit 1; }; \
    done

# hadolint ignore=DL3006
FROM ${MOOSE_MPI_IMAGE}
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG MOOSE_JOBS=4
ARG MARCH=x86-64-v3

# The base image records its environment in /moose_env.sh and
# /opt/mpich/env.sh, which nothing sources outside of apptainer, so set it
# here. MPICH is used rather than OpenMPI, as for MOOSE's own Ubuntu images.
ENV PATH=/opt/mpich/bin:/opt/gcc/bin:/opt/valgrind/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=/opt/mpich/lib \
    MANPATH=/opt/mpich/share/man \
    MOOSE_MPICH_DIR=/opt/mpich \
    CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Packages from apptainer/mpi.def
# hadolint ignore=DL3008
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf automake bison cmake diffutils flex libtool ninja-build patch \
        patchelf patchutils libcppunit-dev libfftw3-dev libglpk-dev libgsl-dev \
        libio-compress-perl libjson-perl libjson-pp-perl libpng-dev libtirpc3 \
        libtirpc-dev zlib1g-dev libc6-dbg && \
    if mpiexec --version | grep hwloc > /dev/null; then \
        apt-get install -y --no-install-recommends hwloc-nox; \
    fi && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /opt/openmpi

# PETSc, following apptainer/petsc.def
# hadolint ignore=DL3003
RUN --mount=type=bind,from=source,source=/moose,target=/root/build/moose,rw \
    cd /root/build/moose && \
    OPTFLAGS="-O2 -march=${MARCH} -mtune=generic" && \
    PETSC_PREFIX=/opt/petsc ./scripts/update_and_rebuild_petsc.sh \
        --skip-submodule-update \
        --with-cc="$CC" --with-cxx="$CXX" --with-fc="$FC" \
        COPTFLAGS="$OPTFLAGS" CXXOPTFLAGS="$OPTFLAGS" FOPTFLAGS="$OPTFLAGS" && \
    cd petsc && \
    make SLEPC_DIR=/opt/petsc PETSC_DIR=/opt/petsc PETSC_ARCH= check
ENV PETSC_DIR=/opt/petsc

# libMesh, following apptainer/libmesh.def (without VTK)
RUN --mount=type=bind,from=source,source=/moose,target=/root/build/moose,rw \
    OPTFLAGS="-march=${MARCH} -mtune=generic" && \
    LIBMESH_DIR=/opt/libmesh METHODS=opt \
        /root/build/moose/scripts/update_and_rebuild_libmesh.sh \
        --skip-submodule-update \
        CFLAGS="$OPTFLAGS" CXXFLAGS="$OPTFLAGS" FCFLAGS="$OPTFLAGS" \
        libmesh_CFLAGS="$OPTFLAGS" libmesh_CXXFLAGS="$OPTFLAGS" libmesh_FCFLAGS="$OPTFLAGS"
ENV LIBMESH_DIR=/opt/libmesh

# WASP, following apptainer/moose-dev.def: built with the underlying compilers
# rather than the MPI wrappers
RUN --mount=type=bind,from=source,source=/moose,target=/root/build/moose,rw \
    CC=gcc CXX=g++ FC=gfortran F90=gfortran F77=gfortran WASP_PREFIX=/opt/wasp \
        /root/build/moose/scripts/update_and_rebuild_wasp.sh \
        --skip-submodule-update \
        -DCMAKE_CXX_FLAGS="-march=${MARCH} -mtune=generic"
ENV WASP_DIR=/opt/wasp \
    PATH=/opt/wasp/bin:${PATH}

# Conduit and MFEM, following apptainer/moose-dev.def
RUN --mount=type=bind,from=source,source=/moose,target=/root/build/moose,rw \
    OPTFLAGS="-march=${MARCH} -mtune=generic" && \
    CONDUIT_DIR=/opt/conduit \
        /root/build/moose/scripts/update_and_rebuild_conduit.sh \
        --skip-submodule-update \
        -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_CXX_FLAGS="$OPTFLAGS" && \
    MPIR_CVAR_ENABLE_GPU=0 CONDUIT_DIR=/opt/conduit MFEM_DIR=/opt/mfem METHODS=opt \
        /root/build/moose/scripts/update_and_rebuild_mfem.sh \
        --skip-submodule-update \
        -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_CXX_FLAGS="$OPTFLAGS"
ENV CONDUIT_DIR=/opt/conduit \
    MFEM_DIR=/opt/mfem

# Python environment, following apptainer/moose-dev.def, at the python version
# that it pins
RUN --mount=type=bind,from=source,source=/moose,target=/root/build/moose \
    PYTHON_VERSION=$(sed -n 's/^ *PYTHON_VERSION=\([0-9.]*\) *$/\1/p' \
        /root/build/moose/apptainer/moose-dev.def) && \
    test -n "$PYTHON_VERSION" && \
    mkdir -p /opt/uv/bin /opt/uv/data && \
    curl -LsSf https://astral.sh/uv/install.sh | \
        env XDG_DATA_HOME=/opt/uv/data XDG_BIN_HOME=/opt/uv/bin UV_NO_MODIFY_PATH=1 sh && \
    export XDG_DATA_HOME=/opt/uv/data && \
    /opt/uv/bin/uv venv /opt/python --no-project --no-config --relocatable --no-cache \
        --python "$PYTHON_VERSION" && \
    UV_PYTHON_DIR="$(realpath "$(dirname "$(realpath "/opt/python/bin/python${PYTHON_VERSION%.*}")")/..")" && \
    ln -s "${UV_PYTHON_DIR}/bin/python3-config" \
        "${UV_PYTHON_DIR}/bin/python${PYTHON_VERSION%.*}-config" /opt/python/bin/ && \
    ln -s "${UV_PYTHON_DIR}/include" /opt/python/ && \
    UV_PYTHON=/opt/python UV_NO_CACHE=1 UV_LINK_MODE=copy \
        /opt/uv/bin/uv pip install -r /root/build/moose/requirements.txt
ENV PATH=/opt/python/bin:${PATH}
