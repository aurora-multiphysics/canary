#!/bin/bash
#SBATCH --job-name=moose_cuda_build_ampere
#SBATCH --output=moose_cuda_build_log.txt
#SBATCH --error=moose_cuda_build_err.txt
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --partition=ampere
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00

## This script builds MOOSE with MFEM configured for CUDA support. It is designed for the Ampere partition on the CSD3 machine
## To run, substitute <ACCOUNT_NAME> with your project account before submitting the job

export METHOD=dbg
export MOOSE_JOBS=32
export STDCXX_PATH=/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-rocky8-x86_64_v3/gcc-8.5.0/gcc-11.3.0-r7adil6umxm6geju4jjmtwo4ecbkzpeo/lib64

module purge
module load rhel8/slurm
module use /usr/local/software/spack/spack-modules/rocky8-a100-20230831/linux-rocky8-zen3
module load openmpi/4.1.5/gcc/33z33ovn

git clone https://github.com/spack/spack.git
. spack/share/spack/setup-env.sh

spack install cmake
spack install python
spack install py-pandas^/$(spack find --format '{hash}' python)
spack install py-pyaml
spack install py-packaging
spack install py-jinja2
spack install py-setuptools
spack install py-deepdiff
spack install py-xmltodict

spack load cmake
spack load python
spack load py-pandas
spack load py-pyaml
spack load py-packaging
spack load py-jinja2
spack load py-setuptools
spack load py-deepdiff
spack load py-xmltodict

git clone https://github.com/Heinrich-BR/moose.git moose_cuda
cd moose_cuda

./scripts/update_and_rebuild_petsc.sh --with-cuda --with-cuda-arch=80
./scripts/update_and_rebuild_libmesh.sh
./scripts/update_and_rebuild_conduit.sh
./scripts/update_and_rebuild_wasp.sh
./scripts/update_and_rebuild_mfem.sh -DMFEM_USE_CUDA=YES -DCUDA_ARCH=sm_80 -DCMAKE_BUILD_TYPE=Debug

./configure --with-mfem

cd framework
cp contrib/mfem/build-dbg/config/config.mk contrib/mfem/build-dbg/config/config.mk.bak
cp contrib/mfem/build-dbg/config/config-install.mk contrib/mfem/build-dbg/config/config-install.mk.bak
cp contrib/mfem/installed/share/mfem/config.mk contrib/mfem/installed/share/mfem/config.mk.bak
sed -i.bak 's/\$<\$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler *>//g' contrib/mfem/build-dbg/config/config.mk
sed -i.bak 's/\$<\$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler *>//g' contrib/mfem/build-dbg/config/config-install.mk
sed -i.bak 's/\$<\$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler *>//g' contrib/mfem/installed/share/mfem/config.mk
make -j $MOOSE_JOBS
cd ../test
make -j $MOOSE_JOBS

export LD_LIBRARY_PATH=${STDCXX_PATH}:$LD_LIBRARY_PATH

