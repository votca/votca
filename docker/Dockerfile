FROM votca/buildenv:latest

ARG TESTOPTS
ARG COVERAGE
ARG CC
ARG CXX
ARG CXXFLAGS
ARG CMAKE_BUILD_TYPE
ARG TESTING
ARG DOXYGEN
ARG MINIMAL
ARG MODULE_BUILD
ARG BUILD_GROMACS
ARG REGRESSION_TESTING

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
RUN rm -rf /home/votca/.ccache
COPY ccache/ /home/votca/.ccache
USER root
RUN chown -R votca:votca /home/votca/votca /home/votca/.ccache
USER votca

WORKDIR /home/votca/votca
RUN mkdir build
WORKDIR build
RUN ccache -z
RUN cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=/usr \
  -DENABLE_TESTING=${TESTING} -DBUILD_CSG_MANUAL=ON -DBUILD_XTP_MANUAL=ON \
  ${COVERAGE:+-DENABLE_COVERAGE_BUILD=ON} -DBUILD_OWN_GROMACS=${BUILD_GROMACS} \
  -DENABLE_REGRESSION_TESTING=${REGRESSION_TESTING} \
  -DMODULE_BUILD=${MODULE_BUILD} ${MODULE_BUILD:+-DCMAKE_INSTALL_PREFIX=/home/votca/votca -DVOTCA_TEST_OPTS="${TESTOPTS}"} \
  ${MINIMAL:+-DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON -DWITH_FFTW=OFF -DWITH_GSL=OFF -DCMAKE_DISABLE_FIND_PACKAGE_GSL=ON \
    -DWITH_MKL=OFF -DCMAKE_DISABLE_FIND_PACKAGE_MKL=ON -DBUILD_MANPAGES=OFF -DWITH_GMX=OFF -DWITH_SQLITE3=OFF \
    -DCMAKE_DISABLE_FIND_PACKAGE_SQLITE3=ON -DBUILD_XTP=OFF -DBUILD_XTP_MANUAL=OFF} \
  ..
RUN make -j2
RUN ccache -s
RUN if [ ${DOXYGEN} ]; then \
      make doxygen; \
      python -m coverxygen --xml-dir share/doc/xml/ --src-dir .. --output doxygen.coverage.info; \
      codecov -R ${PWD}/../${TRAVIS_REPO_SLUG#*/} -F doxygen -X gcov -f doxygen.coverage.info; \
      rm -f doxygen.coverage.info; \
    fi
RUN if [ ${TESTING} = ON ]; then make test CTEST_OUTPUT_ON_FAILURE=1 ARGS="${TESTOPTS}"; fi
RUN if [ ! ${MODULE_BUILD} ]; then \
      make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install/usr && rmdir ${PWD}/install; \
    fi
RUN if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    codecov -R ${PWD}/../${TRAVIS_REPO_SLUG#*/} -F "${CC}" -x "llvm-cov gcov" 2>&1 | grep -v "arcs.*block"; \
  else \
    codecov -R ${PWD}/../${TRAVIS_REPO_SLUG#*/} -F "${CC}" 2>&1 | grep -v "arcs.*block"; \
  fi; \
fi
USER root
RUN make install
USER votca
RUN ccache -s
WORKDIR ../csg-tutorials
