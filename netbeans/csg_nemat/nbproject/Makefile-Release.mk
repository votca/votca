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
FC=g95
AS=

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
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/tools/csg_nemat.o

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
LDLIBSOPTIONS=-L../../../../gmx/install/lib ../../../libcsg/trunk/netbeans/../dist/libcsg.a -lgmx -lboost_program_options -lfftw3 -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/csg_nemat

dist/Release/GNU-Linux-x86/csg_nemat: ../../../libcsg/trunk/netbeans/../dist/libcsg.a

dist/Release/GNU-Linux-x86/csg_nemat: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/csg_nemat ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/tools/csg_nemat.o: nbproject/Makefile-${CND_CONF}.mk ../../src/tools/csg_nemat.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/tools
	${RM} $@.d
	$(COMPILE.cc) -O2 -I../../../libcsg/trunk -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/tools/csg_nemat.o ../../src/tools/csg_nemat.cc

# Subprojects
.build-subprojects:
	cd ../../../libcsg/trunk/netbeans && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/csg_nemat

# Subprojects
.clean-subprojects:
	cd ../../../libcsg/trunk/netbeans && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
