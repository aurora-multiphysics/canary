#!/bin/bash
# =============================================================================
# 02_build_moose_computenode.sh
#
# Build PETSc, libMesh, WASP, Conduit and MFEM (with GSLIB) on a compute
# node, offline. All tarballs and submodules must have been
# pre-fetched by 01_fetch_submodules_loginnode.sh first.
#
# There are two ways to run this script:
# 1: Interactive (request a node first):
#               srun --pty --time=08:00:00 --ntasks=8 --mem=32G bash
#               bash 02_build_moose_computenode.sh
#
# 2:Batch job submission (submit the script to SLURM):
#               sbatch 02_build_moose_computenode.sh
#
# SLURM DIRECTIVES - edit as needed:
#SBATCH --job-name=moose_build
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=moose_build_%j.log
#SBATCH --error=moose_build_%j.err
# =============================================================================

set -e
set -u

# =============================================================================
# CONFIGURATION - edit paths here if your layout differs
# =============================================================================
MOOSE_DIR="$HOME/moose"
PKGS_DIR="$HOME/petsc_packages"  # tarballs downloaded by script 01
PETSC_DIR="$MOOSE_DIR/petsc"
PETSC_ARCH="arch-moose"
GSLIB_DIR="$PKGS_DIR/gslib"     # pre-cloned by script 01
NPROCS="${SLURM_NTASKS:-8}"

# Output helpers
info() { echo "[INFO] $*"; }
warn() { echo "[WARN] $*"; }
error() {
    echo "[ERROR] $*"
    exit 1
}
step() {
    echo ""
    echo "=== $* ==="
    echo "=================================================="
}
ok() { echo "[OK] $*"; }

check() {
    local desc="$1" test_cmd="$2"
    if eval "$test_cmd" >/dev/null 2>&1; then
        ok "$desc"
    else
        echo "[MISSING] $desc"
        ERRORS=$((ERRORS + 1))
    fi
}

# SECTION 0 - Pre-build sanity checks
# Confirm is is running on a compute node and all inputs from script 01 are pointed correctly
step "Pre-flight checks"

info "Running on: $(hostname)"

# Refuse to run on a login node
if echo "$(hostname)" | grep -qi "login"; then
    error "Must run on a COMPUTE node, not a login node.
    Request one: srun --pty --time=08:00:00 --ntasks=8 --mem=32G bash"
fi

# Compute nodes have no internet - warn only if internet unexpectedly present
if curl -s --max-time 3 https://github.com >/dev/null 2>&1; then
    warn "This node has internet access (unusual for Cumulus compute nodes)."
fi

ERRORS=0

# MOOSE repo
check "MOOSE directory ($MOOSE_DIR)"   "[ -d '$MOOSE_DIR' ]"
check "MOOSE configure script"         "[ -f '$MOOSE_DIR/configure' ]"
check "update_and_rebuild_mfem.sh"     "[ -f '$MOOSE_DIR/scripts/update_and_rebuild_mfem.sh' ]"

# PETSc submodule
check "PETSc submodule ($PETSC_DIR)"   "[ -d '$PETSC_DIR' ]"
check "PETSc configure script"         "[ -f '$PETSC_DIR/configure' ]"

# Key tarballs in PKGS_DIR - check by glob pattern
check "HDF5 tarball"         "ls '$PKGS_DIR'/hdf5*.tar.gz 2>/dev/null | head -1 | grep -q ."
check "MUMPS tarball"        "ls '$PKGS_DIR'/MUMPS*.tar.gz 2>/dev/null | head -1 | grep -q ."
check "ScaLAPACK tarball"    "ls '$PKGS_DIR'/scalapack*.tar.gz 2>/dev/null | head -1 | grep -q ."
check "SuperLU_dist tarball" "ls '$PKGS_DIR'/superlu*.tar.gz 2>/dev/null | head -1 | grep -q ."
check "STRUMPACK tarball"    "ls '$PKGS_DIR'/v8.0.0.tar.gz '$PKGS_DIR'/STRUMPACK*.tar.gz 2>/dev/null | head -1 | grep -q ."
# HYPRE tarball: PETSc names it after the commit hash (40 hex chars), not *hypre*
# Accept either a friendly name or the hash-named file (starts with hex digits)
check_hypre() {
    find "$PKGS_DIR" -maxdepth 1 -type f \( \
        -name "*hypre*" -o -name "*HYPRE*" \
    \) 2>/dev/null | head -1 | grep -q . && return 0
    # Hash-named: 40 lowercase hex chars + .tar.gz
    find "$PKGS_DIR" -maxdepth 1 -type f -name "*.tar.gz" 2>/dev/null |
        grep -qE "/[0-9a-f]{40}\.tar\.gz$"
}
check "HYPRE tarball" "check_hypre"

# Submodules that needed manual checkout in script 01
check "libMesh submodule"    "[ -f '$MOOSE_DIR/libmesh/configure.ac' ]"
# Conduit: check for CMakeLists.txt OR src/ directory (both indicate a valid checkout)
check "Conduit checked out"  "[ -f '$MOOSE_DIR/framework/contrib/conduit/CMakeLists.txt' ] || [ -d '$MOOSE_DIR/framework/contrib/conduit/src' ]"
check "MFEM checked out"     "[ -f '$MOOSE_DIR/framework/contrib/mfem/CMakeLists.txt' ]"
check "GSLIB pre-cloned"     "[ -d '$GSLIB_DIR/.git' ]"

[ "$ERRORS" -eq 0 ] || error "$ERRORS pre-flight check(s) failed.
    Run 01_fetch_submodules_loginnode.sh on the login node first."

info "All pre-flight checks passed."

# =============================================================================
# SECTION 1 - Load modules
# Same toolchain used to pre-fetch in script 01
# =============================================================================
step "Section 1: Loading modules"

module purge
module load GCC/13.3.0
module load OpenMPI/5.0.3-GCC-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load flex/2.6.4-GCCcore-13.3.0
module load libtirpc/1.3.5-GCCcore-13.3.0
module load git/2.45.1-GCCcore-13.3.0

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77
export PETSC_DIR="$PETSC_DIR"
export PETSC_ARCH="$PETSC_ARCH"

# MOOSE_JOBS controls parallelism inside the rebuild scripts
export MOOSE_JOBS="$NPROCS"
export METHODS=opt
export METHOD=opt

# Verify all required tools are reachable after module load
for cmd in mpicc mpicxx mpif90 cmake flex git python3; do
    command -v "$cmd" >/dev/null 2>&1 || error "$cmd not found after module load."
done

info "Modules loaded:"
module list 2>&1 | tail -n +2 || true
info "Compiler : $(mpicc --version | head -1)"
info "CMake    : $(cmake --version | head -1)"

# SECTION 2 - PETSc
# Locate HDF5 tarball automatically by glob so the filename does not need
# to be hardcoded here.
# ptscotch, kokkos, kokkos-kernels and umpire are disabled - see build log.
step "Section 2: Building PETSc"

HDF5_TARBALL=$(ls "$PKGS_DIR"/hdf5*.tar.gz 2>/dev/null | head -1 || true)
[ -n "$HDF5_TARBALL" ] || error "HDF5 tarball not found in $PKGS_DIR"
info "HDF5 tarball: $HDF5_TARBALL"

cd "$MOOSE_DIR"
./scripts/update_and_rebuild_petsc.sh \
    --skip-submodule-update \
    --with-packages-download-dir="$PKGS_DIR" \
    --download-ptscotch=0 \
    --download-kokkos=0 \
    --download-kokkos-kernels=0 \
    --download-umpire=0 \
    --download-hdf5="$HDF5_TARBALL"

ok "PETSc"

# SECTION 3 - libMesh
step "Section 3: Building libMesh"

cd "$MOOSE_DIR"
./scripts/update_and_rebuild_libmesh.sh

ok "libMesh"

# SECTION 4 - WASP
step "Section 4: Building WASP"

cd "$MOOSE_DIR"
./scripts/update_and_rebuild_wasp.sh

ok "WASP"

# SECTION 5 - Conduit
step "Section 5: Building Conduit"

cd "$MOOSE_DIR"
./scripts/update_and_rebuild_conduit.sh

ok "Conduit"

# SECTION 6 - MFEM with GSLIB via update_and_rebuild_mfem.sh 9 (tricky part)
#
# GSLIB was pre-cloned by script 01 into $PKGS_DIR/gslib to avoid the
# compute-node git clone failure that happens when MFEM_FETCH_GSLIB=YES
# tries to reach github.com during make.
#
# The script reads GSLIB_SOURCE_DIR to find an existing local clone instead
# of fetching from the internet.
#
# PETSc_CONFIG_CURRENT=YES is exported to skip the multipass PETSc link-test
# that hangs for hours on Lustre filesystems (see build log for details).

step "Section 6: Building MFEM (with GSLIB)"

# Verify PETSc actually built before proceeding
ARCH_DIR="$PETSC_DIR/$PETSC_ARCH"
[ -f "$ARCH_DIR/lib/libpetsc.so" ] \
    || error "PETSc library not found at $ARCH_DIR/lib/libpetsc.so - did Section 2 succeed?"

# Point the MFEM build script at the pre-cloned GSLIB
export GSLIB_SOURCE_DIR="$GSLIB_DIR"

# Skip the multipass PETSc link-test that hangs on Lustre (see build log)
export PETSc_CONFIG_CURRENT=YES

cd "$MOOSE_DIR"
./scripts/update_and_rebuild_mfem.sh

ok "MFEM"

# SECTION 7 - MOOSE configure --with-mfem
step "Section 7: MOOSE configure --with-mfem"

cd "$MOOSE_DIR"
./configure --with-mfem

ok "MOOSE configure"

# SECTION 8 - Build MOOSE modules/combined
step "Section 8: Building MOOSE (modules/combined)"

cd "$MOOSE_DIR/modules/combined"
make METHOD=opt -j"$NPROCS"

ok "MOOSE modules/combined"

# SECTION 9 - Build tests
step "Section 9: Building MOOSE tests"

cd "$MOOSE_DIR/test"
make METHOD=opt -j"$NPROCS"

ok "MOOSE tests"

# FINAL SUMMARY
echo ""
echo "============================================================"
echo "  MOOSE/MFEM build complete on $(hostname)"
echo "  $(date)"
echo "============================================================"
echo ""
echo "  Key install paths:"
echo "    PETSc   : $ARCH_DIR"
echo "    libMesh : $MOOSE_DIR/libmesh/installed"
echo "    MFEM    : $MOOSE_DIR/framework/contrib/mfem/installed"
echo "    MOOSE   : $MOOSE_DIR/modules/combined/combined-opt"
echo ""
echo "  To run the MFEM tests:"
echo "    cd $MOOSE_DIR/test"
echo "    ./run_tests --opt -j$NPROCS --re mfem"
echo ""
echo "  To use MOOSE in a job script, load:"
echo "    module load GCC/13.3.0 OpenMPI/5.0.3-GCC-13.3.0"
echo "    export PETSC_DIR=$PETSC_DIR"
echo "    export PETSC_ARCH=$PETSC_ARCH"
echo "============================================================"