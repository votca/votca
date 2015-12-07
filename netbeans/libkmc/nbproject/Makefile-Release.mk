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
CND_PLATFORM=GNU-Linux
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
	${OBJECTDIR}/_ext/1ce05ac7/bsumtree.o \
	${OBJECTDIR}/_ext/1ce05ac7/event.o \
	${OBJECTDIR}/_ext/1ce05ac7/events.o \
	${OBJECTDIR}/_ext/1ce05ac7/gnode.o \
	${OBJECTDIR}/_ext/1ce05ac7/graphcubic.o \
	${OBJECTDIR}/_ext/1ce05ac7/graphkmc.o \
	${OBJECTDIR}/_ext/1ce05ac7/kmcapplication.o \
	${OBJECTDIR}/_ext/1ce05ac7/kmccalculatorfactory.o \
	${OBJECTDIR}/_ext/1ce05ac7/kmclifetime.o \
	${OBJECTDIR}/_ext/1ce05ac7/longrange.o \
	${OBJECTDIR}/_ext/1ce05ac7/nodedevice.o \
	${OBJECTDIR}/_ext/1ce05ac7/profile.o \
	${OBJECTDIR}/_ext/1ce05ac7/statereservoir.o \
	${OBJECTDIR}/_ext/1ce05ac7/version.o \
	${OBJECTDIR}/_ext/1ce05ac7/version_nb.o \
	${OBJECTDIR}/_ext/1ce05ac7/vssmgroup.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a

${OBJECTDIR}/_ext/1ce05ac7/bsumtree.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/bsumtree.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/bsumtree.o ../../src/libkmc/bsumtree.cc

${OBJECTDIR}/_ext/1ce05ac7/event.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/event.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/event.o ../../src/libkmc/event.cc

${OBJECTDIR}/_ext/1ce05ac7/events.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/events.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/events.o ../../src/libkmc/events.cc

${OBJECTDIR}/_ext/1ce05ac7/gnode.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/gnode.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/gnode.o ../../src/libkmc/gnode.cc

${OBJECTDIR}/_ext/1ce05ac7/graphcubic.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/graphcubic.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/graphcubic.o ../../src/libkmc/graphcubic.cc

${OBJECTDIR}/_ext/1ce05ac7/graphkmc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/graphkmc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/graphkmc.o ../../src/libkmc/graphkmc.cc

${OBJECTDIR}/_ext/1ce05ac7/kmcapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/kmcapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/kmcapplication.o ../../src/libkmc/kmcapplication.cc

${OBJECTDIR}/_ext/1ce05ac7/kmccalculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/kmccalculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/kmccalculatorfactory.o ../../src/libkmc/kmccalculatorfactory.cc

${OBJECTDIR}/_ext/1ce05ac7/kmclifetime.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/kmclifetime.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/kmclifetime.o ../../src/libkmc/kmclifetime.cc

${OBJECTDIR}/_ext/1ce05ac7/longrange.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/longrange.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/longrange.o ../../src/libkmc/longrange.cc

${OBJECTDIR}/_ext/1ce05ac7/nodedevice.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/nodedevice.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/nodedevice.o ../../src/libkmc/nodedevice.cc

${OBJECTDIR}/_ext/1ce05ac7/profile.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/profile.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/profile.o ../../src/libkmc/profile.cc

${OBJECTDIR}/_ext/1ce05ac7/statereservoir.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/statereservoir.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/statereservoir.o ../../src/libkmc/statereservoir.cc

${OBJECTDIR}/_ext/1ce05ac7/version.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/version.o ../../src/libkmc/version.cc

${OBJECTDIR}/_ext/1ce05ac7/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/version_nb.o ../../src/libkmc/version_nb.cc

${OBJECTDIR}/_ext/1ce05ac7/vssmgroup.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libkmc/vssmgroup.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce05ac7
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../include -I../../include -I../../include/votca/kmc -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce05ac7/vssmgroup.o ../../src/libkmc/vssmgroup.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibkmc.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
