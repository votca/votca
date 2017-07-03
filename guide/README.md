[Based on the 2011-2015 CECAM meetings]

Please pick topics you are most interested in, since finishing the tutorial might take longer than one afternoon.
Installation

Building executables

The simplest way is to use the build script which can be downloaded from github project
prefix=~/votca # installation folder (executables will be in the $prefix/bin)
mkdir -p ${prefix}/src # source code folder
cd ${prefix}/src
wget https://raw.githubusercontent.com/votca/buildutil/master/build.sh # get the build script
chmod +x build.sh
./build.sh --prefix ${prefix} -ud tools csg # install tools and csg
The following option might be needed to disable sqlite3 libraries -DWITH_SQLITE3=OFF. 
To build a gromacs version use
./build.sh --prefix ${prefix} -ud gromacs
Downloading the tutorials

git clone https://github.com/votca/csg-tutorials.git csg-tutorials
or            

./build.sh --prefix ${prefix} -du csg-tutorials
Downloading the manual

The corresponding version of the manual can be found on github (development version here), also have a look at the online reference:
http://www.votca.org/Documentation/coarse-graining
Running GROMACS

If you are not familiar with the molecular dynamics (MD) package GROMACS, the following section provides a brief overview on how to modify and run MD simulations using GROMACS, which is used in the tutorial for generating the reference data. You can find all input files for a box of atomistic water in csg-tutorials/spce/atomistic
Input files

You will need four files to run MD simulations
conf.gro - stores the coordinates of the molecule(s). Can be viewed with vmd. The default  file name is conf.gro. 
grompp.mdp stores all simulation options, such as time step, number of simulation steps, etc. 
topol.top - topology of the molecule
forcefield.itp - description of the atomistic force-field (not always needed)
During this tutorial, you will need to modify 
nsteps - number of MD steps
nstxout - output frequency of coordinates to the trajectory file traj.trr
nstfout - output frequency of forces to traj.trr
nstxtcout - output frequency to the traj.xtc file, often used by the iterative Boltzmann inversion method
MD simulations

To run MD simulations using GROMACS, one first must create a binary topology file topol.tpr using the grompp program and then run the MD integrator using the mdrun program. topol.tpr contains conf.gro, grompp.mdp, topol.top, and forcefield.itp 
grompp -c input.gro # -c is needed if other than conf.gro file is used
mdrun -v # -v (verbose) gives the estimate of the run time
Running other MD programs

In addition to GROMACS VOTCA supports ESPResSo, Lammps, dl_poly and ESPResSo++. The interface to these is a bit more advanced, meaning VOTCA will allow you to more crazy things and warn you less about settings, which might not make sense at all. Let have a look at csg-tutorials/spce/ibi_espresso. (ibi_lammps, ibi_espressopp, ibi_dlpoly, are pretty similar)
Input files

You will some files to run MD simulations

spce.gro - stores the coordinates of the molecule(s). Can be viewed with vmd. A pdb or xyz file would be okay, too.
spce.tcl is the simulation script, which will be called by csg_inverse. It stores the whole the simulation procedure.
topol.xml - topology of the molecule, defined in the initial condition. This is needed as most gro/pdb/xyz files have not molecule definition in them. 
Mapping an atomistic onto a coarse-grained trajectory

Water

Folder: csg-tutorials/spce/atomistic/
Have a look at the center of mass mapping file water.xml in this folder. Use csg_map to create a coarse-grained configuration. Visualize both configurations with vmd.
Create a mapping file where the center of a coarse-grained water molecule is the center of charge of he atomistic one.
Hexane

Folder: csg-tutorials/hexane/atomistic/
Run a few steps of atomistic simulations.Take a look at the mapping file, hexane.xml. Calculate the coarse-grained rdf and bond distribution using csg_stat and csg_boltzmann from the atomistic trajectory. csg_stat needs an additional settings file. Use csg-tutorials/hexane/ibi/settings.xml. Visualize the result using vmd.
DPPC (advanced)

Get the files from the plumx project and create a center of geometry/charge mapping. Visualize the result using vmd.
Iterative Boltzmann inversion for SPC/E water

Here an one-site coarse-grained (CG) model of a rigid 3-site water molecule (SPC/E model) is constructed using the iterative Boltzmann inversion (IBI) method. The center of the CG bead is chosen to be the center of mass (COM) of a molecule. The radial distribution function (RDF) is then calculated using the bead coordinates. In the last step a coarse-grained potential is obtained by matching the RDFs of the atomisitc and CG systems using the IBI method.
Atomistic simulations

Run a short MD simulation of a box of SPC/E water using GROMACS. Input files can be found in the folder spce/atomistic of the votca tutorials. Due to limited time, decrease the number of steps (nsteps) in grompp.mdp to a reasonable value (5000)
grompp # combines conf.gro, topol.top, and grompp.mdp, and forcefield.itp into topol.tpr
mdrun -v # runs MD integrator. The trajectory is saved to traj.xtc 
Mapping

Check the mapping file water.xml. Atom names listed in the definition of the COM bead should correspond to those used in the conf.gro file. Use csg_dump to check this
csg_dump --top topol.tpr 
Copy the setting file fmatch.xml to dist.xml and remove all fmatch options (2 blocks) from dist.xml. After this, calculate the center of mass RDF using csg_stat
csg_stat --top topol.tpr --trj traj.trr --cg water.xml --options dist.xml
Compare your RDF of CG-CG.dist.tgt in spce/ibi
Running IBI

Reduce the number of MD steps in grompp.mdp and start the IBI iterations
csg_inverse --options settings.xml
Calculate the pressure after several iterations using g_energy. Apply pressure correction. Add a post update to the settings file, so that a pressure is applied.

If you are done - have a look at the hexane tutorial.

Relative entropy minimization for SPC/E water

Relative entropy (RE) minimization-based coarse-graining of SPC/E water is similar to the IBI example above. The reference atomistic simulations and mapping are the same as in the IBI example.
The water-water CG potential is modeled using a cubic B-spline functional form. An initial guess for the cubic B-spline knot values is provided as CG-CG.param.init. At each iteration step, the CG potential table is computed from the current CG parameters (CG-CG.param.cur), the CG-MD simulation is performed, and CG-CG RDF (CG-CG.dist.new) is determined. Finally, the new CG potential parameters (CG-CG.param.new) are computed using the relative entropy minimization algorithm.
Running RE

Reduce the number of MD steps in grompp.mdp and start the RE iterations
csg_inverse --options settings.xml

          If you are done - have a look at the methanol-water tutorial.

Force matching for liquid methanol

We will now derive the non-bonded CG potentials for liquid methanol using the force matching method. 
Atomistic simulations

Generate the atomistic trajectory. You will need an equilibrated configuration to start from. The rest of the input files can be prepared using the previous tutorial (methanol/atomistic). Make sure that GROMACS outputs both coordinates and forces (nstxout and nstfout should have the same value in grompp.mdp). Due to the limited time of this tutorial, only a short reference run can be sampled. 1000 steps MD run with output every 50 steps should provide reasonable statistics. 
Generate the reference trajectory (traj.trr), type
grompp -c methanol_final.gro
mdrun -v
Running force matching

To prepare the force-matching input, add the force-matching block to the settings file. To estimate min and max use csg_stat with the existing settings.xml to calculate the RDF. 
min and max have to be within the range of the RDF
Change the spline grid, blocksize and parameter constrainedLS. Due to the lack of statistics, the results will be noisy, but you should get a feeling for the whole procedure.
Additional topics

If you are curious you can also compare CG potentials/RDFs obtained using the IBI and force matching methods for SPC/E water. You will see that a single site water model with a pair interaction potential either reproduces the RDF of the atomistic reference, or a chosen thermodynamic property. (Reason: three-body contributions are important but cannot be projected on a coarse-grained force-field).
Hybrid force matching for hexane

In some cases, force matching can have problems with bonded interactions, especially if the functional form of the coarse-grained force field lacks essential interactions such as bond-angle or 3-body correlations. In such cases in can help to perform force matching only on the non-bonded contributions as was shown here.

The files for the tutorial can be found in hexane/hybrid_force_matching. The folder should contain all necessary files to reproduce the plots from the publication. Try to understand the input files and reproduce some of the key figures.

Atomistic simulations

The files to generate the reference run can be found in the folder md. The length of the reference simulation was reduced to a minimum and can be regenerated in a few minutes. The script run.sh also performs a full force matching run for bonded and non-bonded interactions. Compare these results to the hybrid approach.

Force matching for the non-bonded interactions

To be able to only treat non-bonded interactions via force matching, all other contributions need to be excluded. This is achieved by generating a second topology, where all bonded interactions were deleted. Furthermore, all intramolecular interactions were explicitly excluded. The forces of the reference trajectory can then be recalculated using the gromacs mdrun rerun functionality

mdrun -rerun traj.trr

The main modification to the grompp.mdp file is that forces and positions are written out in every interaction, and therefore each configuration from the reference trajectory is written to the output. Then, force-matching can be performed on this newly obtained trajectory.

Regularization of the inverse Monte Carlo method

Inverse Monte Carlo (IMC) needs a well defined cross-correlation matrix for which enough sampling is needed. If there is not enough sampling the algorithm might not converge to a stable solution. This might also happen if the initial potential guess for the iterative scheme is too far away from the real solution of the inverse problem. To overcome this deficiency and to stabilize the algorithm one could apply the so called Tikhonov regularization, which is a common technique to regularize ill-posed inverse problems. For further information on the Tikhonov regularization and/or ill-posed inverse problems in general don't hesitate to have a look at the manual of VOTCA-1.4 to get a short overview or for a more detailed description at this publication or consult any book of choice on regularization of inverse problems. 

This tutorial can be considered to be a proof of concept. It is based on the above mentioned publication. Here the user should get familiar with the application of the Tikhonov regularization and should see its benefit. The file run.sh will execute a preliminary run of 10 steps of iterative Boltzmann inverson (IBI) before the IMC method is applied. The users should figure out what happens if the preliminary IBI steps are skipped and should test different regularization parameters (e.g. 10,100.300,1000). The folder also contains a short python script which performs a singular value decomposition of the cross-correlation matrix (svd.py). Based on this decomposition one could get an educated guess on the order of the magnitude of the regularization parameter. It should be larger than the smallest singular values squared and smaller compared to the larger ones.
Advanced topics

Extending the scripting framework

Write a post update script, which smooths the tail of a potential by transforming dU(r) to s(r)dU(r) with 
s(r) = 1 for r < rt
s(r) = 1-(rc-rt)-3(r - rt)2(3rc-rt-2r) for rt < r < rc
s(r) = 0 for r > rt
Hints:
Start from skeleton.pl and use pressure_cor_simple.pl as a template.
Writing an analysis tool

VOTCA allows to write your own analysis code. There are many examples and two templates for serial and threaded analysis. If you are willing to learn how to write your own analysis in C++, ask for assistance. 

