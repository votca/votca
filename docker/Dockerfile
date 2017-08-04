from votca/buildenv:latest

ARG TESTOPTS
ARG COVERAGE
ARG CC
ARG CXX
ARG CXXFLAGS
ARG CMAKE_BUILD_TYPE

#for coverage
ARG CI
ARG TRAVIS
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT
ARG TRAVIS_OS_NAME

COPY votca/ /home/votca/votca
COPY ccache/ /home/votca/.ccache
USER root
RUN chown -R votca:votca /home/votca/votca /home/votca/.ccache
USER votca

WORKDIR /home/votca/votca
RUN mkdir build
WORKDIR build
RUN ccache -z
RUN cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DBUILD_CTP=ON -DCMAKE_INSTALL_PREFIX=/usr -DLIB=lib64 \
  -DENABLE_TESTING=ON -DBUILD_CSG_MANUAL=ON -DBUILD_XTP_MANUAL=ON  ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} ..
RUN make VERBOSE=1 -j2
RUN make test CTEST_OUTPUT_ON_FAILURE=1 ARGS="${TESTOPTS}"
RUN make install DESTDIR=${PWD} && rm -rf ${PWD}/install
RUN if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-root . --root ../${TRAVIS_REPO_SLUG#*/} -F "${CC}" --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov --gcov-root . --root ../${TRAVIS_REPO_SLUG#*/} -F "${CC}"; \
  fi; \
fi
USER root
RUN make install
USER votca
RUN ccache -s
WORKDIR ../csg-tutorials
