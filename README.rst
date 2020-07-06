|Codacy Badge| |codecov| |CI| |DOI|

This repository contains some small analysis programs which are written
based on the votca csg framework

Programs
========

-  *radii*: calculate gyration radius + hydrodynamic radius of a
   molecule or a set of molecules
-  *sphericalorder*: calculates a spherical order parameter i.e. the
   distribution of e\_r\*u where e\_r is the unit vector from a
   reference molecule to the solvent molecules and u is the principle
   axis of inertia of the solvent molecule
-  *fluctuations*: calculate the number density fluctuations in
   subvolumes. Subvolumes are either cubic slabs along one of the sim
   box axis or speherical subvolumes around a solute
-  *orientcorr*: calculates the distance dependent orientational
   correlation function of a polymer melt
-  *part\_dist*: outputs the time-averaged number of particles, listed
   by particle types (was a part of csg before)
-  *partial\_rdf*: calculates the rdf in a spherical subvolume
-  *traj\_force*: add/subtracts reference forces from a given trajectory
   and stores in a new trajectory

To add your own program just create a new folder and put your \*.cc
files in there

.. |Codacy Badge| image:: https://app.codacy.com/project/badge/Grade/c22e2fb9c1f24dfa9afaed66ecfd7c08
   :target: https://www.codacy.com/gh/votca/csgapps?utm_source=github.com&utm_medium=referral&utm_content=votca/csgapps&utm_campaign=Badge_Grade
.. |codecov| image:: https://codecov.io/gh/votca/csgapps/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/votca/csgapps
.. |CI| image:: https://github.com/votca/csgapps/workflows/CI/badge.svg
   :target: https://github.com/votca/csgapps/actions?query=branch%3Agithub_actions+workflow%3ACI
.. |DOI| image:: https://zenodo.org/badge/23939372.svg
   :target: https://zenodo.org/badge/latestdoi/23939372
