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
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/484457893/apolarsite.o \
	${OBJECTDIR}/_ext/484457893/calculatorfactory.o \
	${OBJECTDIR}/_ext/238600121/eanalyze.o \
	${OBJECTDIR}/_ext/238600121/eimport.o \
	${OBJECTDIR}/_ext/238600121/einternal.o \
	${OBJECTDIR}/_ext/238600121/emultipole_stdal.o \
	${OBJECTDIR}/_ext/238600121/eoutersphere.o \
	${OBJECTDIR}/_ext/238600121/ianalyze.o \
	${OBJECTDIR}/_ext/238600121/idft.o \
	${OBJECTDIR}/_ext/238600121/iimport.o \
	${OBJECTDIR}/_ext/238600121/izindo.o \
	${OBJECTDIR}/_ext/238600121/jobwriter.o \
	${OBJECTDIR}/_ext/238600121/neighborlist.o \
	${OBJECTDIR}/_ext/238600121/pairdump.o \
	${OBJECTDIR}/_ext/238600121/qmmm.o \
	${OBJECTDIR}/_ext/238600121/rates.o \
	${OBJECTDIR}/_ext/238600121/sandbox.o \
	${OBJECTDIR}/_ext/238600121/stateserver.o \
	${OBJECTDIR}/_ext/238600121/tdump.o \
	${OBJECTDIR}/_ext/484457893/ewald2d.o \
	${OBJECTDIR}/_ext/484457893/ewald3d.o \
	${OBJECTDIR}/_ext/484457893/ewaldnd.o \
	${OBJECTDIR}/_ext/484457893/fragment.o \
	${OBJECTDIR}/_ext/484457893/job.o \
	${OBJECTDIR}/_ext/484457893/molecule.o \
	${OBJECTDIR}/_ext/484457893/orbitals.o \
	${OBJECTDIR}/_ext/484457893/parallelpaircalc.o \
	${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o \
	${OBJECTDIR}/_ext/484457893/pewald3d.o \
	${OBJECTDIR}/_ext/484457893/polarseg.o \
	${OBJECTDIR}/_ext/484457893/polarsite.o \
	${OBJECTDIR}/_ext/484457893/polartop.o \
	${OBJECTDIR}/_ext/484457893/progressobserver.o \
	${OBJECTDIR}/_ext/484457893/qmapplication.o \
	${OBJECTDIR}/_ext/484457893/qmcalculator.o \
	${OBJECTDIR}/_ext/484457893/qmdatabase.o \
	${OBJECTDIR}/_ext/484457893/qmmachine.o \
	${OBJECTDIR}/_ext/484457893/qmnblist.o \
	${OBJECTDIR}/_ext/484457893/qmpackagefactory.o \
	${OBJECTDIR}/_ext/484457893/qmpair.o \
	${OBJECTDIR}/_ext/484457893/qmtool.o \
	${OBJECTDIR}/_ext/484457893/segment.o \
	${OBJECTDIR}/_ext/484457893/segmenttype.o \
	${OBJECTDIR}/_ext/484457893/statesaversqlite.o \
	${OBJECTDIR}/_ext/484457893/toolfactory.o \
	${OBJECTDIR}/_ext/484457893/topology.o \
	${OBJECTDIR}/_ext/484457893/version.o \
	${OBJECTDIR}/_ext/484457893/version_nb.o \
	${OBJECTDIR}/_ext/484457893/xinductor.o \
	${OBJECTDIR}/_ext/484457893/xinteractor.o \
	${OBJECTDIR}/_ext/484457893/xjob.o \
	${OBJECTDIR}/_ext/484457893/xmapper.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${OBJECTDIR}/_ext/484457893/apolarsite.o: ../../src/libctp/apolarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/apolarsite.o ../../src/libctp/apolarsite.cc

${OBJECTDIR}/_ext/484457893/calculatorfactory.o: ../../src/libctp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/calculatorfactory.o ../../src/libctp/calculatorfactory.cc

${OBJECTDIR}/_ext/238600121/eanalyze.o: ../../src/libctp/calculators/eanalyze.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eanalyze.o ../../src/libctp/calculators/eanalyze.cc

${OBJECTDIR}/_ext/238600121/eimport.o: ../../src/libctp/calculators/eimport.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eimport.o ../../src/libctp/calculators/eimport.cc

${OBJECTDIR}/_ext/238600121/einternal.o: ../../src/libctp/calculators/einternal.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/einternal.o ../../src/libctp/calculators/einternal.cc

${OBJECTDIR}/_ext/238600121/emultipole_stdal.o: ../../src/libctp/calculators/emultipole_stdal.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/emultipole_stdal.o ../../src/libctp/calculators/emultipole_stdal.cc

${OBJECTDIR}/_ext/238600121/eoutersphere.o: ../../src/libctp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eoutersphere.o ../../src/libctp/calculators/eoutersphere.cc

${OBJECTDIR}/_ext/238600121/ianalyze.o: ../../src/libctp/calculators/ianalyze.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/ianalyze.o ../../src/libctp/calculators/ianalyze.cc

${OBJECTDIR}/_ext/238600121/idft.o: ../../src/libctp/calculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/idft.o ../../src/libctp/calculators/idft.cc

${OBJECTDIR}/_ext/238600121/iimport.o: ../../src/libctp/calculators/iimport.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/iimport.o ../../src/libctp/calculators/iimport.cc

${OBJECTDIR}/_ext/238600121/izindo.o: ../../src/libctp/calculators/izindo.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/izindo.o ../../src/libctp/calculators/izindo.cc

${OBJECTDIR}/_ext/238600121/jobwriter.o: ../../src/libctp/calculators/jobwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/jobwriter.o ../../src/libctp/calculators/jobwriter.cc

${OBJECTDIR}/_ext/238600121/neighborlist.o: ../../src/libctp/calculators/neighborlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/neighborlist.o ../../src/libctp/calculators/neighborlist.cc

${OBJECTDIR}/_ext/238600121/pairdump.o: ../../src/libctp/calculators/pairdump.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/pairdump.o ../../src/libctp/calculators/pairdump.cc

${OBJECTDIR}/_ext/238600121/qmmm.o: ../../src/libctp/calculators/qmmm.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/qmmm.o ../../src/libctp/calculators/qmmm.cc

${OBJECTDIR}/_ext/238600121/rates.o: ../../src/libctp/calculators/rates.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/rates.o ../../src/libctp/calculators/rates.cc

${OBJECTDIR}/_ext/238600121/sandbox.o: ../../src/libctp/calculators/sandbox.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/sandbox.o ../../src/libctp/calculators/sandbox.cc

${OBJECTDIR}/_ext/238600121/stateserver.o: ../../src/libctp/calculators/stateserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/stateserver.o ../../src/libctp/calculators/stateserver.cc

${OBJECTDIR}/_ext/238600121/tdump.o: ../../src/libctp/calculators/tdump.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/tdump.o ../../src/libctp/calculators/tdump.cc

${OBJECTDIR}/_ext/484457893/ewald2d.o: ../../src/libctp/ewald2d.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/ewald2d.o ../../src/libctp/ewald2d.cc

${OBJECTDIR}/_ext/484457893/ewald3d.o: ../../src/libctp/ewald3d.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/ewald3d.o ../../src/libctp/ewald3d.cc

${OBJECTDIR}/_ext/484457893/ewaldnd.o: ../../src/libctp/ewaldnd.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/ewaldnd.o ../../src/libctp/ewaldnd.cc

${OBJECTDIR}/_ext/484457893/fragment.o: ../../src/libctp/fragment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/fragment.o ../../src/libctp/fragment.cc

${OBJECTDIR}/_ext/484457893/job.o: ../../src/libctp/job.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/job.o ../../src/libctp/job.cc

${OBJECTDIR}/_ext/484457893/molecule.o: ../../src/libctp/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/molecule.o ../../src/libctp/molecule.cc

${OBJECTDIR}/_ext/484457893/orbitals.o: ../../src/libctp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/orbitals.o ../../src/libctp/orbitals.cc

${OBJECTDIR}/_ext/484457893/parallelpaircalc.o: ../../src/libctp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/parallelpaircalc.o ../../src/libctp/parallelpaircalc.cc

${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o: ../../src/libctp/parallelxjobcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o ../../src/libctp/parallelxjobcalc.cc

${OBJECTDIR}/_ext/484457893/pewald3d.o: ../../src/libctp/pewald3d.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/pewald3d.o ../../src/libctp/pewald3d.cc

${OBJECTDIR}/_ext/484457893/polarseg.o: ../../src/libctp/polarseg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polarseg.o ../../src/libctp/polarseg.cc

${OBJECTDIR}/_ext/484457893/polarsite.o: ../../src/libctp/polarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polarsite.o ../../src/libctp/polarsite.cc

${OBJECTDIR}/_ext/484457893/polartop.o: ../../src/libctp/polartop.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polartop.o ../../src/libctp/polartop.cc

${OBJECTDIR}/_ext/484457893/progressobserver.o: ../../src/libctp/progressobserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/progressobserver.o ../../src/libctp/progressobserver.cc

${OBJECTDIR}/_ext/484457893/qmapplication.o: ../../src/libctp/qmapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmapplication.o ../../src/libctp/qmapplication.cc

${OBJECTDIR}/_ext/484457893/qmcalculator.o: ../../src/libctp/qmcalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmcalculator.o ../../src/libctp/qmcalculator.cc

${OBJECTDIR}/_ext/484457893/qmdatabase.o: ../../src/libctp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmdatabase.o ../../src/libctp/qmdatabase.cc

${OBJECTDIR}/_ext/484457893/qmmachine.o: ../../src/libctp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmmachine.o ../../src/libctp/qmmachine.cc

${OBJECTDIR}/_ext/484457893/qmnblist.o: ../../src/libctp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmnblist.o ../../src/libctp/qmnblist.cc

${OBJECTDIR}/_ext/484457893/qmpackagefactory.o: ../../src/libctp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpackagefactory.o ../../src/libctp/qmpackagefactory.cc

${OBJECTDIR}/_ext/484457893/qmpair.o: ../../src/libctp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpair.o ../../src/libctp/qmpair.cc

${OBJECTDIR}/_ext/484457893/qmtool.o: ../../src/libctp/qmtool.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmtool.o ../../src/libctp/qmtool.cc

${OBJECTDIR}/_ext/484457893/segment.o: ../../src/libctp/segment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/segment.o ../../src/libctp/segment.cc

${OBJECTDIR}/_ext/484457893/segmenttype.o: ../../src/libctp/segmenttype.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/segmenttype.o ../../src/libctp/segmenttype.cc

${OBJECTDIR}/_ext/484457893/statesaversqlite.o: ../../src/libctp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/statesaversqlite.o ../../src/libctp/statesaversqlite.cc

${OBJECTDIR}/_ext/484457893/toolfactory.o: ../../src/libctp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/toolfactory.o ../../src/libctp/toolfactory.cc

${OBJECTDIR}/_ext/484457893/topology.o: ../../src/libctp/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/topology.o ../../src/libctp/topology.cc

${OBJECTDIR}/_ext/484457893/version.o: ../../src/libctp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version.o ../../src/libctp/version.cc

${OBJECTDIR}/_ext/484457893/version_nb.o: ../../src/libctp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version_nb.o ../../src/libctp/version_nb.cc

${OBJECTDIR}/_ext/484457893/xinductor.o: ../../src/libctp/xinductor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xinductor.o ../../src/libctp/xinductor.cc

${OBJECTDIR}/_ext/484457893/xinteractor.o: ../../src/libctp/xinteractor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xinteractor.o ../../src/libctp/xinteractor.cc

${OBJECTDIR}/_ext/484457893/xjob.o: ../../src/libctp/xjob.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xjob.o ../../src/libctp/xjob.cc

${OBJECTDIR}/_ext/484457893/xmapper.o: ../../src/libctp/xmapper.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xmapper.o ../../src/libctp/xmapper.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
