# syntax=docker/dockerfile:1
#
# canary-moose
#
# MOOSE, configured with MFEM support, built on canary-moose-dev and installed
# to /opt/moose. As Canary has no executable of its own, only the framework
# and moose_test are built.
#
# Required build arguments:
#
#   MOOSE_DEV_IMAGE: The canary-moose-dev image built for the same MOOSE_REVISION
#   MOOSE_REVISION: MOOSE commit to build
#

ARG MOOSE_DEV_IMAGE=unset
# hadolint ignore=DL3006
FROM ${MOOSE_DEV_IMAGE}
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG MOOSE_JOBS=4
ARG MOOSE_REVISION

# The source tree is removed once MOOSE is installed, but every object in the
# syntax dump records its .C file under it, so a workflow that builds
# documentation with this image has to mount a checkout of MOOSE_REVISION
# (which the image sets) at the same place. The directory is left in place as
# the mount point.
# hadolint ignore=DL3003
RUN test -n "$MOOSE_REVISION" && \
    git clone --filter=blob:none --no-checkout https://github.com/idaholab/moose.git /opt/src/moose && \
    cd /opt/src/moose && \
    git checkout --quiet "$MOOSE_REVISION" && \
    ./configure --prefix=/opt/moose --with-mfem="$MFEM_DIR" && \
    cd test && \
    make -j "$MOOSE_JOBS" && \
    make install -j "$MOOSE_JOBS" MOOSE_SKIP_DOCS=1 && \
    find /opt/moose/lib -type f -name "*.so*" -exec strip -S {} + && \
    chmod -R o=u-w,g=u-w /opt/moose && \
    rm -rf /opt/src/moose && \
    mkdir -p /opt/src/moose

ENV PATH=/opt/moose/bin:${PATH} \
    PYTHONPATH=/opt/moose/share/moose/python \
    MOOSE_REVISION=${MOOSE_REVISION}
