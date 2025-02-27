# Build stage with Spack pre-installed and ready to be used

FROM spack/rockylinux8:0.23.1 as builder

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment && \
    set -o noclobber \
    &&  (echo spack: \
    &&   echo '  specs:' \
    &&   echo '  - fftw@3 +openmp ~mpi' \
    &&   echo '  - cuda@12' \
    &&   echo '  - gsl' \
    &&   echo '  - "flex@:2.6.2,2.6.4:"' \
    &&   echo '  - bison' \
    &&   echo '  - python@3.11' \
    &&   echo '  - gperftools' \
    &&   echo '  - intel-oneapi-tbb' \
    &&   echo '  - meson' \
    &&   echo '  - git' \
    &&   echo '  - gmake' \
    &&   echo '  - numactl' \
    &&   echo '  - py-pip' \
    &&   echo '  - py-numpy' \
    &&   echo '  - py-numexpr' \
    &&   echo '  - py-numba ^llvm~clang' \
    &&   echo '  - py-scipy' \
    &&   echo '  - py-matplotlib' \
    &&   echo '  - py-contexttimer' \
    &&   echo '  - py-cffi' \
    &&   echo '  - py-scikit-learn' \
    &&   echo '  - py-asdf' \
    &&   echo '  - py-click' \
    &&   echo '  - py-tqdm' \
    &&   echo '  - py-pybind11' \
    &&   echo '  concretizer:' \
    &&   echo '    unify: true' \
    &&   echo '  config:' \
    &&   echo '    install_tree:' \
    &&   echo '      root: /opt/software' \
    &&   echo '    build_stage:' \
    &&   echo '    - /opt/spack-stage' \
    &&   echo '  view: /opt/views/view' \
    &&   echo '  modules:' \
    &&   echo '    prefix_inspections:' \
    &&   echo '      ./bin:' \
    &&   echo '      - PATH' \
    &&   echo '      ./lib:' \
    &&   echo '      - LD_LIBRARY_PATH' \
    &&   echo '      - LIBRARY_PATH' \
    &&   echo '      ./lib64:' \
    &&   echo '      - LD_LIBRARY_PATH' \
    &&   echo '      - LIBRARY_PATH' \
    &&   echo '      ./include:' \
    &&   echo '      - CPATH' \
    &&   echo '      ./man:' \
    &&   echo '      - MANPATH' \
    &&   echo '      ./share/man:' \
    &&   echo '      - MANPATH' \
    &&   echo '      ./share/aclocal:' \
    &&   echo '      - ACLOCAL_PATH' \
    &&   echo '      ./lib/pkgconfig:' \
    &&   echo '      - PKG_CONFIG_PATH' \
    &&   echo '      ./lib64/pkgconfig:' \
    &&   echo '      - PKG_CONFIG_PATH' \
    &&   echo '      ./share/pkgconfig:' \
    &&   echo '      - PKG_CONFIG_PATH' \
    &&   echo '      ./:' \
    &&   echo '      - CMAKE_PREFIX_PATH' \
    ) > /opt/spack-environment/spack.yaml

# Default to 1 process and $(nproc) jobs
ARG SPACK_INSTALL_PROCESSES=1
ARG SPACK_INSTALL_JOBS

# Install the software with a parallel build,
# and remove unnecessary deps
RUN --mount=type=tmpfs,target=/opt/spack-stage,size=32G \
    --mount=type=cache,target=/opt/spack-cache \
    rm -rf /opt/spack-cache/* && \
    spack mirror add --autopush --unsigned jenkins-cache /opt/spack-cache && \
    # buildcache update is allowed to fail, the cache may be empty!
    spack buildcache update-index jenkins-cache; \
    SPACK_INSTALL_JOBS=${SPACK_INSTALL_JOBS:-$(nproc)} && \
    echo "Building with $SPACK_INSTALL_PROCESSES processes and $SPACK_INSTALL_JOBS jobs" && \
    cd /opt/spack-environment && \
    spack env activate . && \
    for _ in $(seq 1 $SPACK_INSTALL_PROCESSES); do \
        spack install --fail-fast -j $SPACK_INSTALL_JOBS & \
        pids+=($!); \
    done && \
    # catch any errors from the parallel installs
    for pid in ${pids[*]}; do \
        wait $pid ; \
    done && \
    spack gc -y

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh

RUN . /opt/spack-environment/activate.sh && \
    pip install -v \
        git+https://github.com/AleksTk/table-logger.git \
        abacusutils

# Bare OS image to run the installed executables
FROM docker.io/rockylinux:9

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

# TODO: just use spack compilers once we have a build cache
# or remove this entirely if we end up building Abacus in the Spack image
RUN --mount=type=cache,target=/var/cache/dnf \
    dnf update -y \
    && dnf install -y \
        gcc \
        gcc-c++ \
        libatomic

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view


LABEL "app"="abacus"
# N.B. this entrypoint only applies to "docker run" and not to Jenkins
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]
