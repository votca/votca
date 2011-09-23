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
	${OBJECTDIR}/_ext/484457893/qmapplication.o \
	${OBJECTDIR}/_ext/238600121/egaussian.o \
	${OBJECTDIR}/_ext/484457893/qmpair.o \
	${OBJECTDIR}/_ext/238600121/ecorrelation.o \
	${OBJECTDIR}/_ext/484457893/qmdatabase.o \
	${OBJECTDIR}/_ext/484457893/version_nb.o \
	${OBJECTDIR}/_ext/484457893/qmtopology.o \
	${OBJECTDIR}/_ext/238600121/ecoulomb.o \
	${OBJECTDIR}/_ext/484457893/statesaversqlite.o \
	${OBJECTDIR}/_ext/484457893/qmnblist.o \
	${OBJECTDIR}/_ext/484457893/calculatorfactory.o \
	${OBJECTDIR}/_ext/238600121/eoutersphere.o


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
	"${MAKE}"  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/liblibctp.a

dist/Release/GNU-Linux-x86/liblibctp.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/liblibctp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/liblibctp.a

${OBJECTDIR}/_ext/484457893/qmapplication.o: ../../src/libctp/qmapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmapplication.o ../../src/libctp/qmapplication.cc

${OBJECTDIR}/_ext/238600121/egaussian.o: ../../src/libctp/calculators/egaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/egaussian.o ../../src/libctp/calculators/egaussian.cc

${OBJECTDIR}/_ext/484457893/qmpair.o: ../../src/libctp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpair.o ../../src/libctp/qmpair.cc

${OBJECTDIR}/_ext/238600121/ecorrelation.o: ../../src/libctp/calculators/ecorrelation.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/ecorrelation.o ../../src/libctp/calculators/ecorrelation.cc

${OBJECTDIR}/_ext/484457893/qmdatabase.o: ../../src/libctp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmdatabase.o ../../src/libctp/qmdatabase.cc

${OBJECTDIR}/_ext/484457893/version_nb.o: ../../src/libctp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version_nb.o ../../src/libctp/version_nb.cc

${OBJECTDIR}/_ext/484457893/qmtopology.o: ../../src/libctp/qmtopology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmtopology.o ../../src/libctp/qmtopology.cc

${OBJECTDIR}/_ext/238600121/ecoulomb.o: ../../src/libctp/calculators/ecoulomb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/ecoulomb.o ../../src/libctp/calculators/ecoulomb.cc

${OBJECTDIR}/_ext/484457893/statesaversqlite.o: ../../src/libctp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/statesaversqlite.o ../../src/libctp/statesaversqlite.cc

${OBJECTDIR}/_ext/484457893/qmnblist.o: ../../src/libctp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmnblist.o ../../src/libctp/qmnblist.cc

${OBJECTDIR}/_ext/484457893/calculatorfactory.o: ../../src/libctp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/calculatorfactory.o ../../src/libctp/calculatorfactory.cc

${OBJECTDIR}/_ext/238600121/eoutersphere.o: ../../src/libctp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eoutersphere.o ../../src/libctp/calculators/eoutersphere.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/liblibctp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
