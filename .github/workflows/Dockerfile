FROM ghcr.io/votca/buildenv/fedora:latest

COPY votca/ /home/votca/votca
RUN sudo chown -R votca:votca /home/votca/votca

WORKDIR /home/votca/votca
RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=/usr \
  -DENABLE_TESTING=ON -DBUILD_CSG_MANUAL=ON -DENABLE_WERROR=ON \
  -DENABLE_REGRESSION_TESTING=ON \
  -B builddir
RUN cmake --build builddir --parallel 2
RUN cd builddir && ctest --output-on-failure -E regression_
RUN sudo cmake --install builddir
