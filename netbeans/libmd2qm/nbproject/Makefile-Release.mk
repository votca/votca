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
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/76882303/dump_atoms.o \
	${OBJECTDIR}/_ext/76882303/generate_nrgs.o \
	${OBJECTDIR}/_ext/76882303/shufflenrg.o \
	${OBJECTDIR}/_ext/76882303/estatics.o \
	${OBJECTDIR}/_ext/76882303/lambdaout.o \
	${OBJECTDIR}/_ext/1716269789/statesaversqlite.o \
	${OBJECTDIR}/_ext/76882303/polymerrates.o \
	${OBJECTDIR}/_ext/76882303/readxml.o \
	${OBJECTDIR}/_ext/1716269789/qmpair.o \
	${OBJECTDIR}/_ext/1716269789/calculatorfactory.o \
	${OBJECTDIR}/_ext/1716269789/version_nb.o \
	${OBJECTDIR}/_ext/1716269789/qmnblist.o \
	${OBJECTDIR}/_ext/1716269789/qmtopology.o \
	${OBJECTDIR}/_ext/76882303/energycorr.o \
	${OBJECTDIR}/_ext/1716269789/qmapplication.o \
	${OBJECTDIR}/_ext/1716269789/statesaver.o


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
	"${MAKE}"  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/liblibmd2qm.a

dist/Release/GNU-Linux-x86/liblibmd2qm.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/liblibmd2qm.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmd2qm.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/liblibmd2qm.a

${OBJECTDIR}/_ext/76882303/dump_atoms.o: ../../src/libmd2qm/calculators/dump_atoms.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/dump_atoms.o ../../src/libmd2qm/calculators/dump_atoms.cc

${OBJECTDIR}/_ext/76882303/generate_nrgs.o: ../../src/libmd2qm/calculators/generate_nrgs.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/generate_nrgs.o ../../src/libmd2qm/calculators/generate_nrgs.cc

${OBJECTDIR}/_ext/76882303/shufflenrg.o: ../../src/libmd2qm/calculators/shufflenrg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/shufflenrg.o ../../src/libmd2qm/calculators/shufflenrg.cc

${OBJECTDIR}/_ext/76882303/estatics.o: ../../src/libmd2qm/calculators/estatics.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/estatics.o ../../src/libmd2qm/calculators/estatics.cc

${OBJECTDIR}/_ext/76882303/lambdaout.o: ../../src/libmd2qm/calculators/lambdaout.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/lambdaout.o ../../src/libmd2qm/calculators/lambdaout.cc

${OBJECTDIR}/_ext/1716269789/statesaversqlite.o: ../../src/libmd2qm/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/statesaversqlite.o ../../src/libmd2qm/statesaversqlite.cc

${OBJECTDIR}/_ext/76882303/polymerrates.o: ../../src/libmd2qm/calculators/polymerrates.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/polymerrates.o ../../src/libmd2qm/calculators/polymerrates.cc

${OBJECTDIR}/_ext/76882303/readxml.o: ../../src/libmd2qm/calculators/readxml.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/readxml.o ../../src/libmd2qm/calculators/readxml.cc

${OBJECTDIR}/_ext/1716269789/qmpair.o: ../../src/libmd2qm/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/qmpair.o ../../src/libmd2qm/qmpair.cc

${OBJECTDIR}/_ext/1716269789/calculatorfactory.o: ../../src/libmd2qm/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/calculatorfactory.o ../../src/libmd2qm/calculatorfactory.cc

${OBJECTDIR}/_ext/1716269789/version_nb.o: ../../src/libmd2qm/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/version_nb.o ../../src/libmd2qm/version_nb.cc

${OBJECTDIR}/_ext/1716269789/qmnblist.o: ../../src/libmd2qm/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/qmnblist.o ../../src/libmd2qm/qmnblist.cc

${OBJECTDIR}/_ext/1716269789/qmtopology.o: ../../src/libmd2qm/qmtopology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/qmtopology.o ../../src/libmd2qm/qmtopology.cc

${OBJECTDIR}/_ext/76882303/energycorr.o: ../../src/libmd2qm/calculators/energycorr.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/76882303
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/76882303/energycorr.o ../../src/libmd2qm/calculators/energycorr.cc

${OBJECTDIR}/_ext/1716269789/qmapplication.o: ../../src/libmd2qm/qmapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/qmapplication.o ../../src/libmd2qm/qmapplication.cc

${OBJECTDIR}/_ext/1716269789/statesaver.o: ../../src/libmd2qm/statesaver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1716269789
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1716269789/statesaver.o ../../src/libmd2qm/statesaver.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/liblibmd2qm.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
