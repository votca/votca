This folder contains some small analysis programs which are written
based on the votca csg framework, but are not useful for the broad
public.

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

To add your own program just create a new folder and put your `\*.cc`
files in there.
