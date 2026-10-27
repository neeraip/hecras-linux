# HEC-RAS 2D on Linux: the HEC-RAS Linux engines plus this repo's Python tools.
#
#   docker build --platform linux/amd64 -t hecras-linux:6.6 .
#   docker build --platform linux/amd64 --build-arg HECRAS_VERSION=7.0.2 -t hecras-linux:7.0.2 .
#
# The engines are x86_64 builds, so build and run for linux/amd64 (Docker
# Desktop emulates it on Apple silicon, slowly). bin/ and libs/ come from
# git LFS; run `git lfs pull` first if they are small pointer files.
FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 python3-numpy python3-scipy python3-h5py python3-gdal \
        python3-matplotlib python3-shapely zip unzip && \
    rm -rf /var/lib/apt/lists/*

ARG HECRAS_VERSION=6.6
COPY bin/${HECRAS_VERSION}/ /opt/hecras/bin/
COPY libs/ /opt/hecras/libs/
COPY ras_preprocess.py hecras_run.py run_simulation.sh /opt/hecras/
COPY tools/ /opt/hecras/tools/
RUN chmod +x /opt/hecras/bin/* /opt/hecras/run_simulation.sh && \
    test "$(stat -c %s /opt/hecras/bin/RasUnsteady)" -gt 1000000 || \
    (echo "bin/${HECRAS_VERSION}/RasUnsteady is not a real binary - run 'git lfs pull'" && exit 1)

ENV HECRAS_VERSION=${HECRAS_VERSION} \
    PATH=/opt/hecras/bin:$PATH \
    LD_LIBRARY_PATH=/opt/hecras/libs:/opt/hecras/libs/rhel_8:/opt/hecras/libs/mkl

WORKDIR /work
