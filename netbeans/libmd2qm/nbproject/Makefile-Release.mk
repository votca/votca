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
FC=
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
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmtopology.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmnblist.o

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
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/liblibmd2qm.a

dist/Release/GNU-Linux-x86/liblibmd2qm.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/liblibmd2qm.a
	${AR} rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmd2qm.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/liblibmd2qm.a

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmtopology.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libmd2qm/qmtopology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmtopology.o ../../src/libmd2qm/qmtopology.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmnblist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libmd2qm/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libmd2qm/qmnblist.o ../../src/libmd2qm/qmnblist.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/liblibmd2qm.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
