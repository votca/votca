|Codacy Badge| |CI| |Docker| |DOI|

This is VOTCA, which provides the two following subpackages:

-  VOTCA-CSG, a library which provides tools to develop coarse-grained
   potentials from atomistic simulation data
-  VOTCA-XTP, a library providing the DFT+GW-BSE method to calculate 
   electronically excited states in single molecules, in molecular materials 
   with quantum-classical embedding, as well as electron, hole, and exciton 
   dynamics in atomistic MD-trajectories.

Basic installation 
###################
(see performance advice in `Install Guide <share/doc/INSTALL.rst>`__)

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'v2023-rc.2'
    git clone -b ${version} https://github.com/votca/votca.git
    cmake -B builddir -S votca -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix}
    cmake --build builddir --parallel <number of cores>
    cmake --build builddir --target install

Using this code via docker
##########################
::

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

More detailed information:
##########################

1. `Installation <share/doc/INSTALL.rst>`__
2. `Further Information <http://www.votca.org>`__
3. `Developers Guide <share/doc/DEVELOPERS_GUIDE.rst>`__
4. `VOTCA\_LANGUAGE\_GUIDE <share/doc/VOTCA_LANGUAGE_GUIDE.rst>`__
5. `Code of Conduct <share/doc/CODE_OF_CONDUCT.rst>`__

You find the VOTCA-CTP repository
`here <https://gitlab.mpcdf.mpg.de/votca/votca>`__

.. |Codacy Badge| image:: https://app.codacy.com/project/badge/Grade/b5567bfcf2c8411a8057c47fa7126781
   :target: https://www.codacy.com/gh/votca/votca?utm_source=github.com&utm_medium=referral&utm_content=votca/votca&utm_campaign=Badge_Grade
.. |CI| image:: https://github.com/votca/votca/workflows/CI/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ACI+branch%3Amaster
.. |Docker| image:: https://github.com/votca/votca/workflows/Docker/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ADocker+branch%3Amaster
.. |DOI| image:: https://zenodo.org/badge/75022030.svg
   :target: https://zenodo.org/badge/latestdoi/75022030

License:
########

Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)


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
