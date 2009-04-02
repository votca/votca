#
# Gererated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran

# Include project Makefile
include Makefile_nb

# Object Directory
OBJECTDIR=build/profile_release/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgmoleculedef.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgengine.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nematicorder.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/exclusionlist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/connectivity.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologymap.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectorywriter.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbwriter.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectoryreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbtopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/imccoefficients.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/growriter.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/beadlist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/molecule.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nblist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pg
CXXFLAGS=-pg

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../src/libcsg/libcsg.a

../../src/libcsg/libcsg.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libcsg
	${RM} ../../src/libcsg/libcsg.a
	${AR} rv ../../src/libcsg/libcsg.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libcsg/libcsg.a

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgmoleculedef.o: ../../src/libcsg/cgmoleculedef.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgmoleculedef.o ../../src/libcsg/cgmoleculedef.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgengine.o: ../../src/libcsg/cgengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgengine.o ../../src/libcsg/cgengine.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nematicorder.o: ../../src/libcsg/nematicorder.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nematicorder.o ../../src/libcsg/nematicorder.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/exclusionlist.o: ../../src/libcsg/exclusionlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/exclusionlist.o ../../src/libcsg/exclusionlist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/connectivity.o: ../../src/libcsg/connectivity.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/connectivity.o ../../src/libcsg/connectivity.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologymap.o: ../../src/libcsg/topologymap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologymap.o ../../src/libcsg/topologymap.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectorywriter.o: ../../src/libcsg/trajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectorywriter.o ../../src/libcsg/trajectorywriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbwriter.o: ../../src/libcsg/modules/io/pdbwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbwriter.o ../../src/libcsg/modules/io/pdbwriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectoryreader.o: ../../src/libcsg/trajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectoryreader.o ../../src/libcsg/trajectoryreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbtopologyreader.o: ../../src/libcsg/modules/io/pdbtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbtopologyreader.o ../../src/libcsg/modules/io/pdbtopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/imccoefficients.o: ../../src/libcsg/imccoefficients.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/imccoefficients.o ../../src/libcsg/imccoefficients.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/growriter.o: ../../src/libcsg/modules/io/growriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/growriter.o ../../src/libcsg/modules/io/growriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/beadlist.o: ../../src/libcsg/beadlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/beadlist.o ../../src/libcsg/beadlist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o: ../../src/libcsg/modules/io/xmltopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o ../../src/libcsg/modules/io/xmltopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o: ../../src/libcsg/modules/io/gmxtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o ../../src/libcsg/modules/io/gmxtopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o: ../../src/libcsg/modules/io/grotopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o ../../src/libcsg/modules/io/grotopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o: ../../src/libcsg/topologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o ../../src/libcsg/topologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o: ../../src/libcsg/modules/io/gmxtrajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o ../../src/libcsg/modules/io/gmxtrajectorywriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o: ../../src/libcsg/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o ../../src/libcsg/topology.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o: ../../src/libcsg/neighbourlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o ../../src/libcsg/neighbourlist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o: ../../src/libcsg/numberdist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o ../../src/libcsg/numberdist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/molecule.o: ../../src/libcsg/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/molecule.o ../../src/libcsg/molecule.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nblist.o: ../../src/libcsg/nblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nblist.o ../../src/libcsg/nblist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o: ../../src/libcsg/map.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o ../../src/libcsg/map.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o: ../../src/libcsg/modules/io/gmxtrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o ../../src/libcsg/modules/io/gmxtrajectoryreader.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/profile_release
	${RM} ../../src/libcsg/libcsg.a

# Subprojects
.clean-subprojects:
