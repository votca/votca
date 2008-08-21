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
OBJECTDIR=build/Release/GNU-Linux-x86

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
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/moleculeinfo.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptrajectoryreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptopologyreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/configuration.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} dist/Release/GNU-Linux-x86/liblibcsg.a

dist/Release/GNU-Linux-x86/liblibcsg.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/liblibcsg.a
	${AR} rv dist/Release/GNU-Linux-x86/liblibcsg.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/liblibcsg.a

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgmoleculedef.o: ../../src/libcsg/cgmoleculedef.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgmoleculedef.o ../../src/libcsg/cgmoleculedef.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgengine.o: ../../src/libcsg/cgengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/cgengine.o ../../src/libcsg/cgengine.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nematicorder.o: ../../src/libcsg/nematicorder.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/nematicorder.o ../../src/libcsg/nematicorder.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/exclusionlist.o: ../../src/libcsg/exclusionlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/exclusionlist.o ../../src/libcsg/exclusionlist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/connectivity.o: ../../src/libcsg/connectivity.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/connectivity.o ../../src/libcsg/connectivity.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologymap.o: ../../src/libcsg/topologymap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologymap.o ../../src/libcsg/topologymap.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectorywriter.o: ../../src/libcsg/trajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectorywriter.o ../../src/libcsg/trajectorywriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbwriter.o: ../../src/libcsg/modules/io/pdbwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbwriter.o ../../src/libcsg/modules/io/pdbwriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectoryreader.o: ../../src/libcsg/trajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/trajectoryreader.o ../../src/libcsg/trajectoryreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbtopologyreader.o: ../../src/libcsg/modules/io/pdbtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/pdbtopologyreader.o ../../src/libcsg/modules/io/pdbtopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/imccoefficients.o: ../../src/libcsg/imccoefficients.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/imccoefficients.o ../../src/libcsg/imccoefficients.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/growriter.o: ../../src/libcsg/modules/io/growriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/growriter.o ../../src/libcsg/modules/io/growriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/moleculeinfo.o: ../../src/libcsg/moleculeinfo.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/moleculeinfo.o ../../src/libcsg/moleculeinfo.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptrajectoryreader.o: ../../src/libcsg/modules/io/mdptrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptrajectoryreader.o ../../src/libcsg/modules/io/mdptrajectoryreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o: ../../src/libcsg/modules/io/xmltopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/xmltopologyreader.o ../../src/libcsg/modules/io/xmltopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o: ../../src/libcsg/modules/io/gmxtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtopologyreader.o ../../src/libcsg/modules/io/gmxtopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o: ../../src/libcsg/modules/io/grotopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/grotopologyreader.o ../../src/libcsg/modules/io/grotopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o: ../../src/libcsg/topologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topologyreader.o ../../src/libcsg/topologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o: ../../src/libcsg/modules/io/gmxtrajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectorywriter.o ../../src/libcsg/modules/io/gmxtrajectorywriter.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o: ../../src/libcsg/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/topology.o ../../src/libcsg/topology.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o: ../../src/libcsg/neighbourlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/neighbourlist.o ../../src/libcsg/neighbourlist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o: ../../src/libcsg/numberdist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/numberdist.o ../../src/libcsg/numberdist.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptopologyreader.o: ../../src/libcsg/modules/io/mdptopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/mdptopologyreader.o ../../src/libcsg/modules/io/mdptopologyreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o: ../../src/libcsg/map.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/map.o ../../src/libcsg/map.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o: ../../src/libcsg/modules/io/gmxtrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/libcsg/../../src/libcsg/modules/io/gmxtrajectoryreader.o ../../src/libcsg/modules/io/gmxtrajectoryreader.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/configuration.o: ../../src/libcsg/configuration.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/libcsg/configuration.o ../../src/libcsg/configuration.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/liblibcsg.a

# Subprojects
.clean-subprojects:
