#!/bin/bash
# =============================================================================
# 01_fetch_submodules_loginnode.sh
#
# Pre-fetch all submodules and external package tarballs that the
# MOOSE/MFEM build needs. Need to run on the login node (internet access)
#
# Hashes and URLs are read automatically from the checked-out repo,
# so this script works with any MOOSE version.
#
# USAGE   : bash 01_fetch_submodules_loginnode.sh
#           (from any directory — it finds MOOSE via $HOME/moose)
# PS ./01_fetch_submodules_loginnode.sh won't work because it needs sudo rights
# =============================================================================

set -e
set -u

MOOSE_DIR="$HOME/moose"
PKGS_DIR="$HOME/petsc_packages" #directory where 02 point at

# Commands for better readable and pretty output formatting
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

# prechecks
step "Pre-flight checks"

command -v git >/dev/null 2>&1 || error "git not found. Load it: module load git"

[ -d "$MOOSE_DIR" ] || error "MOOSE not found at $MOOSE_DIR"
[ -f "$MOOSE_DIR/.gitmodules" ] || error "Not a git repo: $MOOSE_DIR"

if [ -n "${SLURM_JOB_ID:-}" ]; then
    warn "SLURM_JOB_ID is set - you appear to be on a compute node."
    warn "This script needs internet. Continue anyway? [y/N]"
    read -r ans
    [[ "$ans" =~ ^[Yy]$ ]] || exit 1
fi

info "Checking internet connectivity..."
curl -s --max-time 5 https://github.com >/dev/null 2>&1 \
    || error "Cannot reach github.com - run this on a login node."
info "Internet OK."

mkdir -p "$PKGS_DIR"

# download resume and skips if already exists (function)
download() {
    local url="$1" out="$2"
    if [ -f "$out" ] && [ -s "$out" ] && file "$out" | grep -qE "gzip|bzip2|XZ|Zip|tar"; then
        info "Already present: $(basename "$out") ($(du -sh "$out" | cut -f1))"
        return 0
    fi
    info "Downloading $(basename "$out") ..."
    if command -v wget >/dev/null 2>&1; then
        wget -q --show-progress -c "$url" -O "$out" \
            || { warn "wget failed, retrying with curl..."; curl -L -o "$out" "$url"; }
    else
        curl -L -o "$out" "$url"
    fi
    if ! file "$out" | grep -qE "gzip|bzip2|XZ|Zip|tar"; then
        rm -f "$out"
        error "Download failed or returned non-tarball content for: $url"
    fi
    info "  OK: $(du -sh "$out" | cut -f1)"
}

# get hash of a submodule from git (function)
submodule_hash() {
    # Usage: submodule_hash <moose_dir> <submodule_path>
    git -C "$1" ls-tree HEAD "$2" 2>/dev/null | awk '{print $3}'
}

# clone + checkout a submodule that has update=none function
clone_pinned() {
    local path="$1" url="$2" hash="$3" name="$4"
    if [ -f "$path/CMakeLists.txt" ]; then
        local current
        current=$(git -C "$path" rev-parse HEAD 2>/dev/null || echo "unknown")
        if [ "$current" = "$hash" ]; then
            info "$name already at correct commit $hash - skipping."
            return 0
        else
            warn "$name checked out at $current, expected $hash - re-cloning."
        fi
    fi
    info "Cloning $name @ $hash ..."
    rm -rf "$path"
    git clone "$url" "$path"
    cd "$path"
    git fetch --all --quiet
    git checkout "$hash"
    git submodule update --init --recursive
    info "$name cloned OK"
}

# SECTION 1 — PETSc submodule

step "Section 1: PETSc submodule"

cd "$MOOSE_DIR"
git submodule update --init petsc
info "PETSc submodule OK"

# SECTION 2 — Determine PETSc external package URLs from configure.py
# in the from the moose directory
# Here PETSc configure in "offline dry-run" mode on the login node.
# It exits with code 10 and prints the URLs it needs — then automagically
# download them automatically.

step "Section 2: PETSc external package tarballs (auto-detected URLs)"

# Load python3 (may already be available or need a module)
command -v python3 >/dev/null 2>&1 || error "python3 not found — load it: module load Python"

info "Running PETSc configure in dry-run / download-dir mode to detect URLs..."
info "This will exit with code 10 (normal) and print the list of packages."
echo ""

CONFIGURE_OUT=$(
    cd "$MOOSE_DIR/petsc"
    python3 ./configure \
        --with-packages-download-dir="$PKGS_DIR" \
        --with-debugging=no \
        --with-mpi=1 \
        --with-shared-libraries=1 \
        --with-cxx-dialect=C++17 \
        --with-fortran-bindings=0 \
        --with-sowing=0 \
        --download-fblaslapack=1 \
        --download-hypre=1 \
        --download-metis=1 \
        --download-parmetis=1 \
        --download-mumps=1 \
        --download-scalapack=1 \
        --download-superlu_dist=1 \
        --download-slepc=1 \
        --download-strumpack=1 \
        --download-hdf5=1 \
        --download-libceed=1 \
        --download-openblas=1 \
        --download-hpddm=1 \
        --download-ptscotch=0 \
        --download-kokkos=0 \
        --download-kokkos-kernels=0 \
        --download-umpire=0 \
        2>&1 || true # exit 10 is expected
)

echo "$CONFIGURE_OUT" | grep -v "^Executing\|^stdout\|^$" | head -60 || true
echo ""

info "Parsing and downloading packages listed by PETSc configure..."

# Parse URLs from PETSc output.
# Output format (one package per line):
#   pkgname ['git clone https://...', 'https://.../foo.tar.gz']
#  extract only https:// tarball URLs (not git clone lines).
while IFS= read -r line; do
    # Match lines that start with a word then a space then ['
    if echo "$line" | grep -qE "^[a-zA-Z_][a-zA-Z0-9_]* \['"; then
        pkg=$(echo "$line" | awk '{print $1}')

        # Pull out every https URL ending in .tar.gz or .tgz
        # Use grep -oE (extended, no Perl) for broader compatibility
        urls=$(echo "$line" | grep -oE "https://[^']+\.(tar\.gz|tgz)" || true)

        if [ -z "$urls" ]; then
            warn "No tarball URL found for $pkg — may need manual download"
            continue
        fi

        # Prefer ANL mirror if available, otherwise take the first URL
        url=$(echo "$urls" | grep "cels.anl.gov" | head -1)
        [ -z "$url" ] && url=$(echo "$urls" | head -1)

        # Safety: skip if URL is still empty or doesn't look like a URL
        if [ -z "$url" ] || ! echo "$url" | grep -q "^https://"; then
            warn "Skipping $pkg — could not extract a valid URL from: $line"
            continue
        fi

        # Filename is the last path component of the URL
        fname=$(basename "$url")

        # Safety: skip if filenname is empty or looks like a directory
        if [ -z "$fname" ] || [ -d "$PKGS_DIR/$fname" ]; then
            warn "Skipping $pkg — could not determine output filename from URL: $url"
            continue
        fi

        info "Package: $pkg  ->  $fname"
        download "$url" "$PKGS_DIR/$fname"
    fi
done <<<"$CONFIGURE_OUT"

info "All PETSc tarballs downloaded into $PKGS_DIR"
echo ""
ls -lh "$PKGS_DIR"

# SECTION 3 — libMesh submodules (eigen, netgen)
# netgen on cumulus module goes up to version GCCcore-12.3.0
# Here version >13.0.0 is needed aparently
step "Section 3: libMesh submodules"

cd "$MOOSE_DIR"
git submodule update --init --recursive libmesh
info "libMesh submodules OK"

# SECTION 4 — WASP
step "Section 4: WASP submodule"

cd "$MOOSE_DIR"
git submodule update --init --recursive framework/contrib/wasp 2>/dev/null \
    || info "WASP may be handled internally by update_and_rebuild_wasp.sh — continuing."

# SECTION 5 — Conduit  (update = none in .gitmodules — must clone manually)
step "Section 5: Conduit"

CONDUIT_HASH=$(submodule_hash "$MOOSE_DIR" "framework/contrib/conduit")
[ -n "$CONDUIT_HASH" ] || error "Could not read Conduit hash from git ls-tree."
info "Pinned Conduit commit: $CONDUIT_HASH"

# Resolve relative URL from .gitmodules  (../../LLNL/conduit.git)
CONDUIT_URL="https://github.com/LLNL/conduit.git"

clone_pinned \
    "$MOOSE_DIR/framework/contrib/conduit" \
    "$CONDUIT_URL" \
    "$CONDUIT_HASH" \
    "Conduit"

# SECTION 6 — MFEM  (update = none in .gitmodules — must clone manually)
step "Section 6: MFEM"

MFEM_HASH=$(submodule_hash "$MOOSE_DIR" "framework/contrib/mfem")
[ -n "$MFEM_HASH" ] || error "Could not read MFEM hash from git ls-tree."
info "Pinned MFEM commit: $MFEM_HASH"
MFEM_URL="https://github.com/mfem/mfem.git"

clone_pinned \
    "$MOOSE_DIR/framework/contrib/mfem" \
    "$MFEM_URL" \
    "$MFEM_HASH" \
    "MFEM"

# SECTION 7 — GSLIB  (MFEM fetches this at build time; pre-clone to avoid
#              the compute-node network error during make)
step "Section 7: GSLIB (required by MFEM)"

# Read the tag/commit MFEM's CMakeLists expects
MFEM_CMAKE="$MOOSE_DIR/framework/contrib/mfem/CMakeLists.txt"
GSLIB_TAG=$(grep -A5 "gslib\|GSLIB" "$MFEM_CMAKE" 2>/dev/null \
    | grep -iE "GIT_TAG|VERSION" | head -1 \
    | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' || echo "")

GSLIB_DIR="$PKGS_DIR/gslib"
GSLIB_URL="https://github.com/Nek5000/gslib.git"

if [ -d "$GSLIB_DIR/.git" ]; then
    info "GSLIB already cloned at $GSLIB_DIR"
    cd "$GSLIB_DIR"
    git fetch --all --quiet
else
    info "Cloning GSLIB into $GSLIB_DIR ..."
    git clone "$GSLIB_URL" "$GSLIB_DIR"
    cd "$GSLIB_DIR"
fi

if [ -n "$GSLIB_TAG" ]; then
    info "Checking out GSLIB tag/version: $GSLIB_TAG"
    git checkout "v${GSLIB_TAG}" 2>/dev/null \
        || git checkout "$GSLIB_TAG" 2>/dev/null \
        || warn "Could not checkout tag $GSLIB_TAG — using default branch"
else
    warn "Could not auto-detect GSLIB version from MFEM CMakeLists — using latest main"
    git checkout main 2>/dev/null || git checkout master 2>/dev/null || true
fi
info "GSLIB OK  (stored at $GSLIB_DIR)"

# SECTION 8 — HYPRE hash verification
# (PETSc configure already handled the download above, but we double-check
#  the exact commit hash that MFEM's FindHYPRE.cmake requires, if any)

step "Section 8: Verify HYPRE hash"

HYPRE_FILE=$(find "$PKGS_DIR" -maxdepth 1 -name "*hypre*" -type f 2>/dev/null | head -1)
if [ -n "$HYPRE_FILE" ] && [ -s "$HYPRE_FILE" ]; then
    info "HYPRE tarball present: $(basename "$HYPRE_FILE") ($(du -sh "$HYPRE_FILE" | cut -f1))"
else
    # PETSc configure should have got it; if not, grab it explicitly
    warn "HYPRE tarball not detected — fetching manually..."
    HYPRE_HASH=$(grep -r "hypre" "$MOOSE_DIR/petsc/config/BuildSystem/config/packages/HYPRE.py" \
        2>/dev/null | grep -oE '[0-9a-f]{40}' | head -1 || echo "")
    if [ -n "$HYPRE_HASH" ]; then
        info "HYPRE hash from PETSc config: $HYPRE_HASH"
        download \
            "https://github.com/hypre-space/hypre/archive/${HYPRE_HASH}.tar.gz" \
            "$PKGS_DIR/${HYPRE_HASH}.tar.gz"
    else
        warn "Could not auto-detect HYPRE hash — check $PKGS_DIR manually."
    fi
fi

# FINAL SUMMARY
# Print out summary with relevant info for 02_build_moose_computenode.sh
# Also sanity check

echo ""
echo "============================================================"
echo "  Pre-fetch complete!"
echo "============================================================"
echo ""
echo "  Tarballs staging dir : $PKGS_DIR"
echo "  Tarball count        : $(find "$PKGS_DIR" -maxdepth 1 -name '*.tar.gz' | wc -l)"
du -sh "$PKGS_DIR"
echo ""
echo "  Submodule status:"
git -C "$MOOSE_DIR" submodule status petsc libmesh \
    framework/contrib/conduit framework/contrib/mfem 2>/dev/null || true
echo ""
echo "  GSLIB cached at      : $GSLIB_DIR"
echo ""
echo "  Next step:"
echo "    Request a compute node:"
echo "      srun --pty --time=08:00:00 --ntasks=8 --mem=32G bash"
echo "    Then run:"
echo "      bash 02_build_moose_computenode.sh"
echo "============================================================"
