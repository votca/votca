# Votca Tutorial Guide

(Based on the 2011-2017 CECAM meetings)

Please pick topics you are most interested in, since finishing the tutorial might take longer than one afternoon.

# Installation

## Building VOTCA

The simplest way is to use CMake:

```bash
prefix=~/votca
version=master # or 'stable' or 'v1.4.1'
git clone -b ${version} --recursive https://github.com/votca/votca.git
cd votca
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${prefix} ..
make -j5
make install
```

The following option might be needed to disable sqlite3 libraries `-DWITH_SQLITE3=OFF`. 

To build a gromacs version for VOTCA use

```bash
cmake -DBUILD_OWN_GROMACS=ON ..
make -j5
make install
```

## Using the tutorials

All the tutorials are in the installation folder, i.e. `${prefix}/share/votca/csg-tutorials`.

## Downloading the manual

The corresponding version of the manual can be found on [github](https://github.com/votca/csg-manual/releases) (development version [here](http://doc.votca.org/manual.pdf)).

# Running GROMACS

The tutorial uses GROMACS, a molecular dynamics (MD) package, for generating the reference data. If you are not familiar with GROMACS, please read the following section for a brief overview on how to set up and run MD simulations using GROMACS. You can find all the GROMACS input files in the `atomistic` folder of each tutorial, e.g., for a box of atomistic water see `csg-tutorials/spce/atomistic/`.

## Input files

You will need four files to run MD simulations
  * `conf.gro` - stores the coordinates of the molecule(s). It can be viewed with vmd. The default  file name is `conf.gro`. 
  * `grompp.mdp` - stores all simulation options, such as time step, number of simulation steps, etc. 
  * `topol.top` - topology of the molecule
  * `forcefield.itp` - description of the atomistic force-field (not always needed)

During this tutorial, you will need to modify 
  * `nsteps` - number of MD steps
  * `nstxout` - output frequency of coordinates to the trajectory file traj.trr
  * `nstfout` - output frequency of forces to traj.trr
  * `nstxout-compressed` - output frequency to the traj.xtc file, often used by the iterative Boltzmann inversion method
  * `nstlog` - output frequency of the md.log file
  * `nstenergy` - output frequency to the ener.edr file containing all thermodynamic information

## MD simulations

To run MD simulations using GROMACS, first one must create a binary topology file `topol.tpr` using the `gmx grompp` program. Then run the MD integrator using the `gmx mdrun` program. `topol.tpr` contains `conf.gro`, `grompp.mdp`, `topol.top`, and `forcefield.itp`. 

```bash
gmx grompp -c conf.gro # -c is needed if other than conf.gro file is used
gmx mdrun -v # -v (verbose) gives the estimate of the run time
```

# Running other MD programs

In addition to GROMACS, VOTCA supports `ESPResSo`, `Lammps`, `dl_poly`, `HOOMD-blue`, and `ESPResSo++`. The interface to these is a bit more advanced, meaning VOTCA will allow you to do more crazy things and warn you less about settings, which might not make sense at all. Let's have a look at `csg-tutorials/spce/espressopp/`. (`ibi_lammps`, `ibi_espresso`, `ibi_dlpoly`, and `ibi_hoomd-blue` are pretty similar)

Input files

Here is an example of the ESPResSo++ input files to run the CG MD simulations of the Iterative Boltzmann inversion procedure.

  * `conf.gro` - stores the coordinates of the molecule(s). It can be viewed with vmd. A pdb or xyz file would be okay, too.
  * `spce.py` is the simulation script which will be called by csg_inverse. It stores the whole simulation procedure.
  * `topol.xml` -topology of the molecule, defined in the initial condition. This is needed as most gro/pdb/xyz files have no molecule definition in them.


# Mapping an atomistic trajectory onto a coarse-grained trajectory

## SPC/E Water

Folder: `csg-tutorials/spce/atomistic/`

Have a look at the center of mass mapping file `water.xml` in this folder. Use `csg_map` to create a coarse-grained configuration (adjust the corresponding command in the `run.sh` script to map `conf.gro` instead of `confout.gro`). Visualize both configurations with vmd. (execute `vmd conf.gro &`)

Create a mapping file where the center of a coarse-grained water molecule is not the center of mass, but a different one, e.g., the charge center of the atomistic one.

Run a short MD simulation of a box of SPC/E water using GROMACS. (This is needed later for the different coarse-graining (CG) methods.) The input files are in the same folder (`csg-tutorials/spce/atomistic/`). Due to limited time, decrease the number of steps (`nsteps`) in `grompp.mdp` to a reasonable value (`5000-10000`) and adjust the output frequency of the trajectory, .log, and .edr files to, e.g., 50-100 (`nstxout`, `nstfout`, `nstlog`, nstenergy):

```bash
gmx grompp # combines conf.gro, topol.top, and grompp.mdp, and forcefield.itp (which is in the gromos43a1.ff folder of GROMACS) into topol.tpr
gmx mdrun -v # runs MD integrator. The trajectory is saved to traj.trr
```

Map this short atomistic trajectory to a CG trajectory using `csg_map` (adjust the corresponding line in the `run.sh` script to map `traj.trr` onto `traj_cg.trr` instead of `conf(out).gro` to `conf_cg.gro`.)


# Iterative Boltzmann inversion (IBI) for SPC/E water

Here, a one-site coarse-grained (CG) model of a rigid 3-site water molecule (SPC/E model) is constructed (see the previous section of the tutorial) using the iterative Boltzmann inversion (IBI) method. The center of the CG bead is chosen to be the center of mass (COM) of a molecule. The target radial distribution function (RDF) is calculated from the CG bead coordinates obtained by mapping the reference atomistic trajectory. In the last step, a coarse-grained potential is obtained by matching the RDFs of the atomisitc and CG systems using the IBI method. For a more detailed description, look at the following [publication](http://pubs.acs.org/doi/abs/10.1021/ct900369w).

## Atomistic simulation

A short atomistic MD simulation has been already run in the last part of the previous section (see folder:  `csg-tutorials/spce/atomistic/`). You can extract all thermodynamic information (total energy, kinetic energy, pressure, etc.) from the binary `ener.edr` file using `gmx energy`. Running the `Extract_Energies.sh` script does this for you (see the corresponding command line in the `run.sh` script). It creates an additional subfolder `energies` with all thermodynamic information:

```bash
./Extract_Energies.sh $equi #  The argument $equi is optional. If provided, analysis will start at the corresponding time frame (in GROMACS units, ps) (e.g. 1 to 5)
```

## Calculation of RDF

Once again, check the mapping file `water.xml`. Atom names listed in the definition of the COM bead should correspond to those used in the `conf.gro` file. You can use `csg_dump` to check this:

```bash
csg_dump --top topol.tpr 
```

Check the options file `settings.xml`. It contains the section of the corresponding non-bonded interaction (`CG-CG`). After this, calculate the center of mass RDF:

```bash
csg_stat --top topol.tpr --trj traj.trr --cg water.xml --options settings.xml (--nt 3 --begin $equi) # ( ) denotes additional options: --nt # number of threads to run calculation in parallel with more then one threads (e.g. 3), --begin # time frame in GROMACS units (ps) to start analysis (e.g. 1 to 5)
```
Compare your RDF with `CG-CG.dist.tgt` in `csg-tutorials/spce/ibi/` which has been calculated with a much longer atomistic simulation run.


## Running IBI

Now switch to the folder: `csg-tutorials/spce/ibi/`. Reduce the number of MD steps in `grompp.mdp` and adjust the equilibration time in the settings.xml file (`cg.gromacs.equi_time`) to a lower value (time frame in GROMACS units (ps) at which the analysis of the CG trajectory is started in each IBI step). Finally, start the IBI iterations:

```bash
csg_inverse --options settings.xml
```

At each iteration step, the current CG potential is `CG-CG.pot.cur`. Then, the CG-MD simulation is performed and the CG-CG RDF (`CG-CG.dist.new`) is determined. Finally, the CG potential is updated: `CG-CG.pot.new`.

Calculate the pressure after several iterations using `gmx energy`. You can do this by copying the `Extract_Energies.sh` script from the main folder (`csg-tutorials/spce/ibi/`) to the appropriate step folder (`step_xxx`) and executing it (type `./Extract_Energies.sh`). Again, it will create a subfolder energies. You may notice that the pressure of the CG simulation is significantly too high. This can be adjusted by applying a pressure correction. To do so, add an appropriate post update option to the settings file, so that a (simple) pressure correction is applied. You can check the corresponding section of the manual.


# Inverse Monte Carlo (IMC) for SPC/E water

Developing a CG potential with the inverse Monte Carlo (IMC) method works in a similar way as in the IBI example. The IMC procedure, again, requires the coarse-grained RDF of the single bead mapping based on the atomistic simulation as input. Therefore, in `csg-tutorials/spce/imc`, again, the RDF of a long atomistic simulation run is proveded: `CG-CG.dist.tgt`.

It is required to reduce the number of MD steps in `grompp.mdp` and to adjust the equilibration time in the `settings.xml` file (`cg.gromacs.equi_time`) to a lower value. IMC converges faster than IBI, but needs a better statistical sampling in each CG step. Therefore, the number of MD steps should be larger than in the IBI tutorial. Start the IMC iterations:

```bash
csg_inverse --options settings.xml
```

Again, the current CG potential is `CG-CG.pot.cur`, the CG-CG RDF is `CG-CG.dist.new` and the updated CG potential is: `CG-CG.pot.new`.

Calculate the pressure after several iterations using `gmx energy` and compare it with the IBI tutorial. You can do this by copying the `Extract_Energies.sh` script from the main folder (`csg-tutorials/spce/ibi/`) to the appropriate step folder (`step_xxx`) and executing it.


# Relative entropy (RE) minimization for SPC/E water

Relative entropy (RE) minimization based coarse-graining of SPC/E water works similar to the IBI and the IMC example above. The reference atomistic simulation and the CG mapping are the same as in the IBI example. Again, in `csg-tutorials/spce/re/`, the RDF of a long atomistic simulation run is provided: `CG-CG.dist.tgt`. In this tutorial, the water-water CG potential is modeled using a cubic B-spline functional form. An initial guess for the cubic B-spline knot values is provided in CG-CG.param.init. At each iteration step, the CG potential table is computed from the current CG parameters (`CG-CG.param.cur`), the CG-MD simulation is performed, and the CG-CG RDF (`CG-CG.dist.new`) is determined. Finally, the new CG potential parameters (`CG-CG.param.new`) are computed using the relative entropy minimization algorithm. Reduce the number of MD steps in `grompp.mdp`, adjust the equilibration time in the `settings.xml` file (`cg.gromacs.equi_time`) and start the RE iterations:

```bash
csg_inverse --options settings.xml
```

Again, the current CG potential is `CG-CG.pot.cur`, the CG-CG RDF is `CG-CG.dist.new` and the updated CG potential is: `CG-CG.pot.new`.

Calculate the pressure after several iterations using `gmx energy` and compare it with the IBI tutorial. Again, you can do this by copying the `Extract_Energies.sh` script from the main folder (`csg-tutorials/spce/ibi/`) to the appropriate step folder (`step_xxx`) and executing it. For a more detailed description of the RE method, look at the following [publication](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0131754).


# Force matching for SPC/E water

We will now derive a non-bonded CG potential for SPC/E using the force matching method. 

## Atomistic simulation

Basis for the force matching procedure is an atomistic MD simulation. All files are found in the atomistic folder (`csg-tutorials/spce/atomistic/`). If you have done the above tutorials, you have already generated the files of the atomistic md run. If not, then do so. Adjust the number of time steps to a reasonable value (`5000-10000`) and also choose an appropriate output frequency of the trajectory. Make sure, both, coordinates and forces are written to the trajectory file (`nstxout` and `nstfout` should have the same value in `grompp.mdp`).

## Force matching (FM)

All files for running the actual force matching calculation can be found in  `csg-tutorials/spce/force_matching/`. Have a look at the settings file (`fmatch.xml`). In the general force matching section (`cg.fmatch`), the number of frames to read in simultaneously (`frames_per_block`) and the type of LS solver (constrainedLS) are fixed (it is preferred to use constrained LS). The fmatch block of the interaction (`cg.non-bonded.fmatch`) contains the interaction range (`min` and `max`), the step size for the internal spline representation (`step`) and the output step (`out_step`). `min` and `max` have to be adjusted to be within the range of the RDF (see the calculation of RDF section of the IBI tutorial). Run the FM calculation (see also the corresponding line of the `run.sh` script): 

```bash
csg_fmatch --top ../atomistic/topol.tpr --trj ../atomistic/traj.trr --options fmatch.xml --cg water.xml (--begin $equi ) # ( ) denotes additional option: , --begin # time frame in GROMACS units (ps) to start analysis (e.g. 1 to 5)
```

To obtain the CG potential, the CG force has to be integrated. (see the appropriate lines in the `run.sh` script):

```bash
csg_call table integrate CG-CG.force CG-CG.pot # integrates the table
csg_call table linearop CG-CG.pot CG-CG.pot -1 0 # multiplication of all table values by -1 (potential)
```

Change the spline grid (`step`), blocksize, and parameter constrainedLS. his should provide an overview of the whole procedure.

## Running of CG simulation

To run a CG simulation with GROMACS, the potential has to be converted to a potential table, GROMACS can read (`table_CG_CG.xvg`). (Check the inverse section in the `fmatch.xml` for the corresponding options):

```bash
csg_call --options fmatch.xml --ia-name CG-CG --ia-type non-bonded convert_potential gromacs --clean input.pot table_CG_CG.xvg # calls convert_potential gromacs. Unsampled regions for distances smaller than the min value are extrapolated.
```

To run a CG simulation, you will need the `conf.gro`, `topol.top`, `index.ndx` and `grompp.mdp` files. You can use the ones of the ibi tutorial, and adjust the number of timesteps and output settings. Then run the simulation. Afterwards, you can calculate the RDF and thermodynamic data as explained in the IBI tutorial. You can also use the `Extract_Energies.sh` script of the IBI tutorial. When calculating the RDF from the CG simulation, you don't need a mapping file and the --cg option can be omitted.

Compare the CG potential, the RDF and thermodynamics with the ones of the IBI, IMC or RE method (or any other method) and with the atomistic simulation. You will see that different methods lead to significantly different interaction potentials and a single site water model with a pair interaction potential is not capable of reproducing the RDF and thermodynamics at the same time. (Reason: three-body contributions are important but cannot be projected on a 2-body coarse-grained force-field. The incorporation of non-bonded 3-body interactions is work in progress. An extension to analytic non-bonded 3-body interactions will be released soon).


# Visualization of IBI updates

Go to the folder `csg-tutorials/spce/realtime`. Execute the `run.sh` script. 


# Coarse-graining of liquid methanol

In the folder `csg-tutorials/methanol/`, you will find all relevant files to run an atomistic simulation of liquid methanol and obtain CG potentials with the IBI, IMC and FM method. Look at the SPC/E water tutorial to learn how to do this. You can compare the differences of the CG potentials, RDFs and thermodynamics between the different CG schemes and the atomistic simulation to the differences of the SPC/E water simulations. You will see that in the case of methanol, a pair potential is a better approximation to an ideal CG potential as in the SPC/E water case. The reason is that non-bonded 3-body effects are less important. 


# Coarse-graining of liquid hexane

Go to the folder `csg-tutorials/hexane/`. So far, we only considered single bead mappings. Hexane is a small alkane molecule. In this tutorial, a 3 bead CG mapping with one bond type and one angle type is chosen

## Atomistic simulation

Go into the `csg-tutorials/hexane/atomistic/` folder. Have a look at the mapping file `hexane.xml`. The hexane molecule is mapped to 3 beads with two different bead types with two bonds (of the same type) and one angle. You will find all relevant GROMACS input files in the folder. Have a look at the `run.sh` script. Again, adjust the number of time steps and the output frequencies in grompp.mdp and run an atomistic simulation. Extract the thermodynamic information (`./Extract_Energies.sh`) and calculate the 3 different RDFs (A-A, B-B, A-B) and the bond and angle distributions with `csg_stat`. In addition, you can compute the bond and angle distributions with `csg_boltzmann` (see `run_boltzmann.sh`) Compare the distributions to those in `csg-tutorials/hexane/ibi_all` obtained by a significantly longer atomistic MD run. You can map the (final) .gro file of the atomistic simuation to the CG one to get all necessary information for running the IBI procedure.

## IBI for all interactions

Go to the folder `csg-tutorials/hexane/ibi_all`. The folder contains target RDFs and bond and angle distributions from a longer atomistic MD run. Have a look at the `settings.xml` file. It contains the sections for the non-bonded and bonded interactions. The three non-bonded interactions are updated every 3rd iteration step (first  A-A, then B-B and then A-B, etc., see the `cg.non-bonded.inverse.do_potential` section). The bonded interactions are updated every iteration step. Adjust the number of time steps and output frequency in the `grompp.mdp` file and the equilibration time in the `settings.xml` file (`cg.inverse.gromacs.equi_time`) and start the IBI iterations.
Calculate the pressure after several iterations using `gmx energy`. You can do this by copying the `Extract_Energies.sh` script from the main folder (`csg-tutorials/spce/ibi/`) to the appropriate step folder (`step_xxx`) and executing it.

## IBI for non-bonded interactions only

Go to the folder `csg-tutorials/hexane/ibi_nonbonded`. The folder contains the same target RDFs as the `csg-tutorials/hexane/ibi_all` folder. In addition, it contains the tabulated bond and angle potentials (`table_b1.xvg` and `table_a1.xvg`). They are obtained by (non-iterative) Boltzmann inversion of the bond and angle target distribution functions of a longer MD simulation run. You can compute them with `csg_boltzmann` (see the `run_boltzmann.sh` script). A lack of statistics will become most apparent at lower values of the angle. In this case you might want to adjust the min value in the `boltzmann_cmds` file.  Compare your results with the pre-computed ones (`table_b1.xvg` and `table_a1.xvg`) and with the ones in the subfolder `step_001` in `csg-tutorials/hexane/ibi_all`  Again, adjust the number of time steps and equilibration time and start the iteration process.
Calculate the pressure and compare the thermodynamic properties as well as the obtained non-bonded potentials with the ones of the `csg-tutorials/hexane/ibi_all tutorial`. Calculate the bond and angle distributions (This can be done with `csg_stat`, using the `settings.xml` file and the `hexane_cg.xml` mapping file of the `csg-tutorials/hexane/ibi_all` folder.)

## FM for all interactions together

Go to the folder `csg-tutorials/hexane/force_matching`. The folder contains the hexane mapping file with bond and angle interactions (`hexane.xml`) and the force matching options file (`fmatch.xml`). Have a look at both files and the `run.sh` script and start the force matching procedure. Basis is the atomistic trajectory with force output in `csg-tutorials/hexane/atomistic`. Integrate the force output to obtain the potentials and convert them to GROMACS tables. (see the `run.sh` script). Compare the obtained potentials to the IBI potentials. You can run CG simulations with the CG potentials, again, using the `conf.gro`, `topol.top`, `index.ndx` and `grompp.mdp` files from the `csg-tutorials/hexane/ibi_all` folder. Calculate the RDFs, bond and angle distributions.
When comparing the results to the IBI potentials, you will see that in some cases, force matching can have problems with bonded interactions, especially if the functional form of the coarse-grained force field lacks essential interactions such as bond-angle or 3-body correlations. In such cases it can help to perform force matching only on the non-bonded contributions as was shown [here](http://dx.doi.org/10.1002/mats.201100011).

## FM for non-bonded interactions only

The files for the tutorial can be found in `csg-tutorials/hexane/hybrid_force_matching/`. The folder should contain all necessary files to reproduce the plots from the [publication](http://dx.doi.org/10.1002/mats.201100011). To be able to parametrize only the non-bonded interactions via force matching, an atomistic trajectory has to be generated containing only forces contributing to the non-bonded interactions, meaning all other contributions need to be excluded. This is achieved by generating a second atomistic topoly file. Have a look at the topol.top file and compare it to the one in `csg-tutorials/hexane/atomistic/`. All bonded interactions have been deleted. Furthermore, all intramolecular interactions have been explicitly excluded. Generate the binary GROMACS topology file using this `topol.top` file and the `conf.gro` and `grommp.mdp` file of the reference atomistic trajectory (`csg-tutorials/hexane/atomistic/`). Then, generate the trajectory file with excluded bonded interactions using `gmx mdrun` with the `-rerun` option. (Have a look at the corresponding lines of the `run.sh` script). 
Have a look at the `fmatch.xml` file. It now only contains the non-bonded interactions. Start the FM calculation. Afterwards, intergrate the force output and convert the potentials to GROMACS tables (see the `run.sh` script). You can run the CG simulation, using the `conf.gro`, `topol.top`, `index.ndx` and `grompp.mdp` files and the bond and angle potential (`table_b1.xvg` and `table_a1.xvg`) from the `csg-tutorials/hexane/ibi_nonbonded/` folder. Calculate the RDFs, bond and angle distributions and compare the results to the IBI results and FM of all interactions together.


# Regularization of the inverse Monte Carlo method

For this tutorial go to the folder `csg-tutorials/LJ1-LJ2/imc`. Inverse Monte Carlo (IMC) needs a well defined cross-correlation matrix for which enough sampling is needed. If there is not enough sampling the algorithm might not converge to a stable solution. This might also happen if the initial potential guess for the iterative scheme is too far away from the real solution of the inverse problem. To overcome this deficiency and to stabilize the algorithm one could apply the so called Tikhonov regularization, which is a common technique to regularize ill-posed inverse problems. For further information on the Tikhonov regularization and/or ill-posed inverse problems in general don't hesitate to have a look at the manual of VOTCA-1.4 to get a short overview or for a more detailed description at this [publication](http://epjst.epj.org/articles/epjst/abs/2016/07/epjst160120/epjst160120.html) or consult any book of choice on regularization of inverse problems. 

This tutorial can be considered to be a proof of concept. It is based on the above mentioned [publication](http://epjst.epj.org/articles/epjst/abs/2016/07/epjst160120/epjst160120.html). Here the user should get familiar with the application of the Tikhonov regularization and should see its benefit. The file `run.sh` will execute a preliminary run of 10 steps of iterative Boltzmann inverson (IBI) before the IMC method is applied. The users should figure out what happens if the preliminary IBI steps are skipped and should test different regularization parameters (e.g. 10,100.300,1000). The folder also contains a short python script which performs a singular value decomposition of the cross-correlation matrix (`svd.py`). Based on this decomposition one could get an educated guess on the order of the magnitude of the regularization parameter. It should be larger than the smallest singular values squared and smaller compared to the larger ones.


# Additional tutorials

Have a look in the folder `csg-tutorials`. It contains additional tutorials on propane, methanol-water and urea-water mixtures. To do the tutorials, have a look at the corresponding `run.sh` scripts.



# Advanced topics

## Extending the scripting framework

Write a post update script, which smooths the tail of a potential by transforming `dU(r)` to `s(r)dU(r)` with 

```
s(r) = 1 for r < rt
s(r) = 1-(rc-rt)-3(r - rt)2(3rc-rt-2r) for rt < r < rc
s(r) = 0 for r > rt
```

Hints:
Start from `skeleton.pl` and use `pressure_cor_simple.pl` as a template.

## Writing an analysis tool

VOTCA allows to write your own analysis code. There are many examples and two templates for serial and threaded analysis. If you are willing to learn how to write your own analysis in C++, ask for assistance. 

