|Codacy Badge| |CI| |Docker| |FreeBSD| |DOI| |JOSS|

VOTCA (Versatile Object-oriented Toolkit for Coarse-graining Applications)
is a C++ package for multiscale simulations of soft matter:

-  **VOTCA-XTP**: excited states of molecules and molecular materials with
   DFT and many-body Green's functions (GW-BSE), also in quantum-classical
   (QM/MM) embedding with polarizable environments; electronic couplings,
   and charge and exciton transport (kinetic Monte Carlo, time of flight,
   graph random walks) in morphologies from molecular dynamics.
-  **VOTCA-CSG**: coarse-graining from atomistic simulations (force matching,
   relative entropy, inverse Monte Carlo) and analysis of MD trajectories
   from GROMACS, LAMMPS and other packages.

Documentation, manuals and tutorials: https://www.votca.org

Installation
############

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'v2026'
    git clone -b ${version} https://github.com/votca/votca.git
    cmake -B builddir -S votca -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix}
    cmake --build builddir --parallel <number of cores>
    ctest --test-dir builddir # to test the build
    cmake --build builddir --target install
    source ${prefix}/bin/VOTCARC.bash

Requirements, distribution packages and performance advice (MKL, Apple
Accelerate) are in the `installation guide <share/sphinx/INSTALL.rst>`__.
VOTCA is also available as a Docker image::

    docker run -it votca/votca /bin/bash

Getting started
###############

The programs ``xtp_tools``, ``xtp_map``, ``xtp_run`` and ``xtp_parallel``
(XTP) and ``csg_*`` (CSG) print their options with ``--help``; for the XTP
calculators use ``-l`` (list) and ``-d <name>`` (options). The XTP tutorials
(Jupyter notebooks in ``xtp-tutorials``) cover single-molecule GW-BSE, QM/MM
workflows for GROMACS and LAMMPS morphologies, and charge transport.

Citation
########

If you use VOTCA, please cite

  B. Baumeier et al., VOTCA: multiscale frameworks for quantum and
  classical simulations in soft matter, J. Open Source Softw. 9, 6864
  (2024), https://doi.org/10.21105/joss.06864

The XTP programs print the papers of the methods used at the end of each
run; ``--cite <name>`` lists them for a calculator or tool without running it.

Contributing
############

Bugs and questions: `GitHub issues <https://github.com/votca/votca/issues>`__
and the `VOTCA forum <https://groups.google.com/g/votca>`__. See the
`developers guide <share/sphinx/DEVELOPERS_GUIDE.rst>`__, the
`language guide <share/sphinx/VOTCA_LANGUAGE_GUIDE.rst>`__ and the
`code of conduct <share/sphinx/CODE_OF_CONDUCT.rst>`__. The legacy
VOTCA-CTP code is `here <https://gitlab.mpcdf.mpg.de/votca/votca>`__.

.. |Codacy Badge| image:: https://app.codacy.com/project/badge/Grade/b5567bfcf2c8411a8057c47fa7126781
   :target: https://www.codacy.com/gh/votca/votca?utm_source=github.com&utm_medium=referral&utm_content=votca/votca&utm_campaign=Badge_Grade
.. |CI| image:: https://github.com/votca/votca/actions/workflows/continuous-integration-workflow.yml/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ACI+branch%3Amaster
.. |Docker| image:: https://github.com/votca/votca/actions/workflows/docker-build.yml/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ADocker+branch%3Amaster
.. |FreeBSD| image:: https://github.com/votca/votca/actions/workflows/freebsd.yml/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3AFreeBSD+branch%3Amaster
.. |DOI| image:: https://zenodo.org/badge/75022030.svg
   :target: https://zenodo.org/badge/latestdoi/75022030
.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.06864/status.svg
   :target: https://doi.org/10.21105/joss.06864

License
#######

Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)


Licensed under the Apache License, Version 2.0 (the "License"); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

::

       http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
