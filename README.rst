|Codacy Badge| |CI| |Docker| |DOI|

This is VOTCA's next generation build system for CSG and XTP. It allows
you to easily install:

-  VOTCA-CSG, a library which provides tools to develop coarse-grained
   potentials from atomistic simulation data
-  VOTCA-XTP, a library designed to determine electronic properties of
   organic materials from atomistic MD-trajectories.

Usage:

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v1.6.2'
    git clone -b ${version} --recursive https://github.com/votca/votca.git
    cmake -B builddir -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} votca
    cmake --build builddir --parallel <number of cores>
    cmake --build buildir --target install

Using this code via docker:

::

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

For further details see:

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
