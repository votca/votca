[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f9468cb3863b4a83b75eda6250755a3d)](https://www.codacy.com/manual/votca-package/csgapps?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=votca/csgapps&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/votca/csgapps/branch/master/graph/badge.svg)](https://codecov.io/gh/votca/csgapps)
[![Build Status](https://travis-ci.org/votca/csgapps.svg?branch=master)](https://travis-ci.org/votca/csgapps)
[![pipeline status](https://gitlab.com/votca/csgapps/badges/master/pipeline.svg)](https://gitlab.com/votca/csgapps/commits/master)

This repository contains some small analysis programs which are written based on the votca_csg framework

# Programs

*   __radii__: calculate gyration radius + hydrodynamic radius of a molecule or a set of molecules
*   __sphericalorder__: calculates a spherical order parameter i.e. the distribution of e_r\*u where e_r is the unit vector from a reference molecule to the solvent molecules and u is the principle axis of inertia of the solvent molecule
*   __fluctuations__: calculate the number density fluctuations in subvolumes. Subvolumes are either cubic slabs along one of the sim box axis or speherical subvolumes around a solute
*   __orientcorr__: calculates the distance dependent orientational correlation function of a polymer melt
*   __part_dist__: outputs the time-averaged number of particles, listed by particle types (was a part of csg before)
*   __partial_rdf__: calculates the rdf in a spherical subvolume
*   __traj_force__: add/subtracts reference forces from a given trajectory and stores in a new trajectory

To add your own program just create a new folder and put your *.cc files in there
