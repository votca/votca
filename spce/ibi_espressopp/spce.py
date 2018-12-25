#!/usr/bin/env python

# DON'T RUN THIS MANUALLY
# RUN "csg_inverse --options settings.xml" instead

import sys
import espresso
import mpi4py.MPI as MPI

from espresso.tools import decomp

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

x, y, z, Lx, Ly, Lz = espresso.tools.convert.gromacs.read('conf.gro')
num_particles = len(x)

print 'number of particles: ', num_particles

######################################################################

# system
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(box, nodeGrid, rc, skin)

system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

######################################################################

# adding particles
props = ['id', 'type', 'pos', 'v', 'mass']
new_particles = []
mass = 1.0
type = 0
for pid in range(num_particles):
  pos  = espresso.Real3D(x[pid],y[pid],z[pid])
  vel  = espresso.Real3D(0,0,0)
  part = [pid, type, pos, vel, mass]
  new_particles.append(part)
  
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

##########################################################################################

# interaction
vl = espresso.VerletList(system, cutoff=rc)
tabP = espresso.interaction.Tabulated(itype=1, filename='CG_CG.tab', cutoff=rc)
tabI = espresso.interaction.VerletListTabulated(vl)
tabI.setPotential(type1=0, type2=0, potential=tabP)
system.addInteraction(tabI)

##########################################################################################

# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

# thermostat
lT = espresso.integrator.LangevinThermostat(system)
lT.gamma = 5.0
lT.temperature = 2.5
integrator.addExtension(lT)

##########################################################################################

print "equilibrating ..."
for step in range(100):
  integrator.run(100)
  espresso.tools.info(system, integrator)

print "runing ..."

espresso.tools.info(system, integrator)
for step in range(900):
  integrator.run(100)
  espresso.tools.info(system, integrator)
  print 'writing .xyz trajectory...'
  # we are simulating in nm, but xyz files are in angstroms! -> scale= 10.0
  espresso.tools.DumpConfigurations.fastwritexyz('traj.xyz', system, unfolded = False, append = True, velocities = False, scale = 10.0)
  
print "finished"
