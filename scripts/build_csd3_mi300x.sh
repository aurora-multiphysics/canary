#!/bin/bash
#SBATCH --job-name=moose_hip_build_ampere
#SBATCH --output=moose_hip_build_log.txt
#SBATCH --error=moose_hip_build_err.txt
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --partition=ukaea-mi300x
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00

## This script builds MOOSE with MFEM configured for HIP support. It is designed for the MI300X partition on the CSD3 machine
## To run, substitute <ACCOUNT_NAME> with your project account before submitting the job

export METHOD=opt
export MOOSE_JOBS=32
export HIP_ARCH=gfx942

module purge
module load rhel9/default-amdgpu
module load openmpi/5.0.9/llvm-amdgpu-7.1.1/awm5kxiz
module load gcc/15.2.0

export CC=mpicc
export CXX=mpicxx
export FC=mpifort
export HIPCXXFLAGS="-std=c++17"
export HIPFLAGS="-std=c++17"
export NO_STOP_MESSAGE=1
export ROCM_PATH=/usr/local/software/spack/csd3-rhel9/externals/rocm-7.1.1
export LIBRARY_PATH="${ROCM_PATH}/llvm/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${ROCM_PATH}/llvm/lib:${LD_LIBRARY_PATH}"
export PATH="${ROCM_PATH}/llvm/bin:${PATH}"

BUILD_DIR=$(pwd)/moose_hip
export BUILD_DIR

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}" || exit

git clone --depth=1 --branch v1.1.1 https://github.com/spack/spack.git
# shellcheck source=/dev/null
source spack/share/spack/setup-env.sh
spack install gcc@15.2.0
spack env create amd-env
spack env activate amd-env

cat >>"$(spack location -e amd-env)/spack.yaml" <<'EOF'
  packages:
    openmpi:
      externals:
        - spec: openmpi@5.0.9
          modules:
            - openmpi/5.0.9
      buildable: false
EOF

spack add cmake \
	python \
	py-pandas \
	py-pyaml \
	py-packaging \
	py-jinja2 \
	py-setuptools \
	py-deepdiff \
	py-xmltodict \
	openblas \
	metis \
	parmetis \
	libtirpc \
	hdf5+mpi+hl

spack concretize -f
spack install -j ${MOOSE_JOBS}

git clone https://github.com/xiaoyeli/superlu_dist.git
cd superlu_dist || exit
mkdir build
cd build || exit

cmake .. \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=../install \
	-DTPL_ENABLE_HIPLIB=TRUE \
	-DCMAKE_HIP_ARCHITECTURES=${HIP_ARCH} \
	-DTPL_ENABLE_INTERNAL_BLASLIB=OFF \
	-DTPL_BLAS_LIBRARIES="$(spack location -i openblas)/lib/libopenblas.so" \
	-DTPL_ENABLE_LAPACKLIB=ON \
	-DTPL_LAPACK_LIBRARIES="$(spack location -i openblas)/lib/libopenblas.so" \
	-DTPL_PARMETIS_INCLUDE_DIRS="$(spack location -i parmetis)/include;$(spack location -i metis)/include" \
	-DTPL_PARMETIS_LIBRARIES="$(spack location -i parmetis)/lib/libparmetis.so;$(spack location -i metis)/lib/libmetis.so" \
	-DBUILD_SHARED_LIBS=ON

make -j ${MOOSE_JOBS}
make install
SUPERLU_DIR=$(realpath ../install)
export SUPERLU_DIR

cd "${BUILD_DIR}" || exit

git clone https://github.com/hypre-space/hypre.git
cd hypre || exit
mkdir build
cd build || exit

cmake ../src/ \
	-DHYPRE_ENABLE_HIP=ON \
	-DCMAKE_HIP_ARCHITECTURES=${HIP_ARCH} \
	-DHYPRE_ENABLE_GPU_AWARE_MPI=ON \
	-DMPI_INCLUDE_DIR="$(spack location -i openmpi)/include" \
	-DHYPRE_BUILD_UMPIRE=ON \
	-DCMAKE_INSTALL_PREFIX=../install

make -j ${MOOSE_JOBS}
make install
HYPRE_DIR=$(realpath ../install)
export HYPRE_DIR

cd "${BUILD_DIR}" || exit
git clone https://github.com/Heinrich-BR/moose.git
cd moose || exit
git switch HenriqueBR/lld_llvm_fix # Only here until fix gets merged

./scripts/update_and_rebuild_petsc.sh --with-hip \
	--with-hip-arch=${HIP_ARCH} \
	--with-metis-include="$(spack location -i metis)/include" \
	--with-metis-lib="$(spack location -i metis)/lib/libmetis.so" \
	--with-parmetis-include="$(spack location -i parmetis)/include" \
	--with-parmetis-lib="$(spack location -i parmetis)/lib/libparmetis.so" \
	--with-hdf5-include="$(spack location -i hdf5)/include" \
	--with-hdf5-lib="$(spack location -i hdf5)/lib/libhdf5.so" \
	--download-openblas-make-options='DYNAMIC_ARCH=1 FC=gfortran' \
	--download-hdf5=0 \
	--download-metis=0 \
	--download-parmetis=0 \
	--download-kokkos=0 \
	--download-kokkos-kernels=0 \
	--download-strumpack=0 \
	--download-superlu_dist=0 \
	--download-hypre=0 \
	--download-openblas=1 \
	--with-hypre=1 \
	--with-hypre-dir="${HYPRE_DIR}" \
	--with-superlu_dist=1 \
	--with-superlu_dist-dir="${SUPERLU_DIR}" \
	--with-64-bit-indices=0 \
	--CXXFLAGS+="-include cstring" \
	LIBS="-lhiprtc"

LIBTIRPC_DIR=$(spack location -i libtirpc)
export LIBTIRPC_DIR
export CPPFLAGS="-I${LIBTIRPC_DIR}/include/tirpc ${CPPFLAGS}"
export LDFLAGS="-L${LIBTIRPC_DIR}/lib -ltirpc ${LDFLAGS}"
export PKG_CONFIG_PATH="${LIBTIRPC_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH}"

./scripts/update_and_rebuild_libmesh.sh CC=mpicc CXX=mpicxx FC=mpifort F77=mpifort
./scripts/update_and_rebuild_conduit.sh
./scripts/update_and_rebuild_wasp.sh
./scripts/update_and_rebuild_mfem.sh -DMFEM_USE_HIP=YES -DHIP_ARCH=${HIP_ARCH} -DSuperLUDist_DIR="${SUPERLU_DIR}" -DHDF5_DIR="$(spack location -i hdf5)"

./configure --with-mfem

cd framework || exit
make -j ${MOOSE_JOBS}
cd ../test || exit
make -j ${MOOSE_JOBS}
