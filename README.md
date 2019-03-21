[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b946fb9d8a774893a731917cd07cd0b5)](https://www.codacy.com/app/JoshuaSBrown/csgapps?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=votca/csgapps&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/votca/csgapps/branch/master/graph/badge.svg)](https://codecov.io/gh/votca/csgapps)

This repository contains some small analysis programs which are written based on the votca_csg framework

# Programs

## radii

calculate gyration radius + hydrodynamik radius of a molecule or a set of molecules

## sphericalorder

calculates a spherical order parameter i.e. the distribution of e_r*u where e_r is the unit vector from
a reference molecule to the solvent molecules and u is the principle axis of inertia of the solvent molecule

## fluctuations

calculate the number density fluctuations in subvolumes. Subvolumes are either cubic slabs along one of the
sim box axis or speherical subvolumes around a solute

## orientcorr

calculates the distance dependent orientational correlation function of a polymer melt

## part_dist:

outputs the time-averaged number of particles, listed by particle types (was a part of csg before)

## partial_rdf:

calculates the rdf in a spherical subvolume 

## traj_force:

add/subtracts reference forces from a given trajectory and stores in a new trajectory



To add your own program just create a new folder and put your *.cc files in there
