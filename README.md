[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f9468cb3863b4a83b75eda6250755a3d)](https://www.codacy.com/manual/votca-package/csgapps?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=votca/csgapps&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/votca/csgapps/branch/master/graph/badge.svg)](https://codecov.io/gh/votca/csgapps)
[![Build Status](https://travis-ci.org/votca/csgapps.svg?branch=master)](https://travis-ci.org/votca/csgapps)
[![pipeline status](https://gitlab.com/votca/csgapps/badges/master/pipeline.svg)](https://gitlab.com/votca/csgapps/commits/master)
[![DOI](https://zenodo.org/badge/23939372.svg)](https://zenodo.org/badge/latestdoi/23939372)

This repository contains some small analysis programs which are written based on the votca csg framework

# Programs

*   _radii_: calculate gyration radius + hydrodynamic radius of a molecule or a set of molecules
*   _sphericalorder_: calculates a spherical order parameter i.e. the distribution of e\_r\*u where e\_r is the unit vector from a reference molecule to the solvent molecules and u is the principle axis of inertia of the solvent molecule
*   _fluctuations_: calculate the number density fluctuations in subvolumes. Subvolumes are either cubic slabs along one of the sim box axis or speherical subvolumes around a solute
*   _orientcorr_: calculates the distance dependent orientational correlation function of a polymer melt
*   _part_dist_: outputs the time-averaged number of particles, listed by particle types (was a part of csg before)
*   _partial\_rdf_: calculates the rdf in a spherical subvolume
*   _traj\_force_: add/subtracts reference forces from a given trajectory and stores in a new trajectory

To add your own program just create a new folder and put your *.cc files in there
