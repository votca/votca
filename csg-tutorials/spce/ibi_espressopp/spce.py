#!/usr/bin/env python3

# DON'T RUN THIS MANUALLY
# RUN "csg_inverse --options settings.xml" instead

import sys
import espressopp
import mpi4py.MPI as MPI
from espressopp import Real3D
from espressopp.tools import gromacs
from espressopp.tools import decomp
#from espressopp.tools import timers

'''
2180 Particles
Box 4.031^3
langevin thermostat T=2.5 gamma=5
dt=0.002
cutoff=0.9
'''

rc   = 0.9
skin = 0.3
timestep = 0.002

grofile = "conf.gro"
topfile = "" #"topol.top"

# this calls the gromacs parser for processing the top file (and included files) and the conf file
# The variables at the beginning defaults, types, etc... can be found by calling
# gromacs.read(grofile,topfile) without return values. It then prints out the variables to be unpacked
defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, resname, resid, Lx, Ly, Lz =gromacs.read(grofile,topfile)
num_particles = len(x)

print('number of particles: ', num_particles)

######################################################################

# system
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,box, rc, skin)
cellGrid = decomp.cellGrid(box, nodeGrid, rc, skin)

system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

######################################################################

# adding particles
props = ['id', 'pos', 'v', 'type', 'mass']#, 'q']
new_particles = []
mass = 1.0
type = 0
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
           Real3D(.0,.0,.0), type, mass]
           #Real3D(vx[pid], vy[pid], vz[pid]), types[pid], masses[pid], charges[pid]]
    new_particles.append(part)

system.storage.addParticles(new_particles, *props)
system.storage.decompose()

##########################################################################################

# interaction
vl = espressopp.VerletList(system, cutoff=rc)
tabP = espressopp.interaction.Tabulated(itype=1, filename='CG_CG.tab', cutoff=rc)
tabI = espressopp.interaction.VerletListTabulated(vl)
tabI.setPotential(type1=0, type2=0, potential=tabP)
system.addInteraction(tabI)

##########################################################################################

# integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# thermostat
lT = espressopp.integrator.LangevinThermostat(system)
lT.gamma = 5.0
lT.temperature = 2.5
integrator.addExtension(lT)

# set system properties
int_steps = 900
eq_steps = 100
steps_per_int = 100

##########################################################################################

print("equilibrating ...")
espressopp.tools.analyse.info(system, integrator)
for step in range(eq_steps):
  integrator.run(steps_per_int)
  espressopp.tools.analyse.info(system, integrator)

print("runing ...")

espressopp.tools.analyse.info(system, integrator)
for step in range(int_steps):
  integrator.run(steps_per_int)
  espressopp.tools.analyse.info(system, integrator)
  print('writing .xyz trajectory...')
  # we are simulating in nm, but xyz files are in angstroms! -> scale= 10.0
  espressopp.tools.fastwritexyz('traj.xyz', system, velocities = False, unfolded = False, append = True, scale = 10.0)

print("finished")
