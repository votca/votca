from votca/buildenv:latest

COPY votca/ /home/votca/votca
COPY ccache/ /home/votca/.ccache
USER root
RUN chown -R votca:votca /home/votca/votca /home/votca/.ccache
USER votca

WORKDIR /home/votca/votca
RUN mkdir build
WORKDIR build
RUN ccache -z
RUN cmake -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DBUILD_CTP=ON -DCMAKE_INSTALL_PREFIX=/usr -DLIB=lib64 -DENABLE_TESTING=ON ..
RUN make VERBOSE=1 -j2
RUN make install DESTDIR=${PWD}
USER root
RUN make install
USER votca
RUN ccache -s
WORKDIR ../csg-tutorials
