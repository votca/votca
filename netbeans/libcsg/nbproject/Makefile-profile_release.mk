#
# Generated Makefile - do not edit!
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
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=profile_release
CND_DISTDIR=dist

# Include project Makefile
include Makefile_nb

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgmoleculedef.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgengine.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nematicorder.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/exclusionlist.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/connectivity.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologymap.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectorywriter.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbwriter.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectoryreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbtopologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/growriter.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/beadlist.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/xmltopologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtopologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/grotopologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/esptopologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologyreader.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectorywriter.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/imcio.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topology.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/molecule.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblist.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/map.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblistgrid.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/version_nb.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/csgapplication.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectoryreader.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pg
CXXFLAGS=-pg

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-profile_release.mk ../../src/libcsg/libcsg.a

../../src/libcsg/libcsg.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libcsg
	${RM} ../../src/libcsg/libcsg.a
	${AR} rv ../../src/libcsg/libcsg.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libcsg/libcsg.a

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgmoleculedef.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/cgmoleculedef.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgmoleculedef.o ../../src/libcsg/cgmoleculedef.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgengine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/cgengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/cgengine.o ../../src/libcsg/cgengine.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nematicorder.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/nematicorder.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nematicorder.o ../../src/libcsg/nematicorder.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/exclusionlist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/exclusionlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/exclusionlist.o ../../src/libcsg/exclusionlist.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/connectivity.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/connectivity.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/connectivity.o ../../src/libcsg/connectivity.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologymap.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/topologymap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologymap.o ../../src/libcsg/topologymap.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectorywriter.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/trajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectorywriter.o ../../src/libcsg/trajectorywriter.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbwriter.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/pdbwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbwriter.o ../../src/libcsg/modules/io/pdbwriter.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectoryreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/trajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/trajectoryreader.o ../../src/libcsg/trajectoryreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbtopologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/pdbtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/pdbtopologyreader.o ../../src/libcsg/modules/io/pdbtopologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/growriter.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/growriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/growriter.o ../../src/libcsg/modules/io/growriter.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/beadlist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/beadlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/beadlist.o ../../src/libcsg/beadlist.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/xmltopologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/xmltopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/xmltopologyreader.o ../../src/libcsg/modules/io/xmltopologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtopologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/gmxtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtopologyreader.o ../../src/libcsg/modules/io/gmxtopologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/grotopologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/grotopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/grotopologyreader.o ../../src/libcsg/modules/io/grotopologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/esptopologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/esptopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/esptopologyreader.o ../../src/libcsg/modules/io/esptopologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologyreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/topologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topologyreader.o ../../src/libcsg/topologyreader.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectorywriter.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/gmxtrajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectorywriter.o ../../src/libcsg/modules/io/gmxtrajectorywriter.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/imcio.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/imcio.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/imcio.o ../../src/libcsg/imcio.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topology.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/topology.o ../../src/libcsg/topology.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/molecule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/molecule.o ../../src/libcsg/molecule.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/nblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblist.o ../../src/libcsg/nblist.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/map.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/map.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/map.o ../../src/libcsg/map.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblistgrid.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/nblistgrid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/nblistgrid.o ../../src/libcsg/nblistgrid.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/version_nb.o ../../src/libcsg/version_nb.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/csgapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/csgapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/csgapplication.o ../../src/libcsg/csgapplication.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectoryreader.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libcsg/modules/io/gmxtrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -I../../../include -I/usr/include/libxml2 -I../../../../../ruehle/gmx/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libcsg/modules/io/gmxtrajectoryreader.o ../../src/libcsg/modules/io/gmxtrajectoryreader.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/profile_release
	${RM} ../../src/libcsg/libcsg.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
