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
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/csg_resample/../../src/tools/csg_resample.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lboost_program_options -lgsl -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../src/tools/csg_resample

../../src/tools/csg_resample: ${BUILD_SUBPROJECTS}

../../src/tools/csg_resample: ${OBJECTFILES}
	${MKDIR} -p ../../src/tools
	${LINK.cc} -o ../../src/tools/csg_resample ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/csg_resample/../../src/tools/csg_resample.o: ../../src/tools/csg_resample.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/csg_resample/../../src/tools
	$(COMPILE.cc) -g -I../../../include -I../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/netbeans/csg_resample/../../src/tools/csg_resample.o ../../src/tools/csg_resample.cc

# Subprojects
.build-subprojects:
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../src/tools/csg_resample

# Subprojects
.clean-subprojects:
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Release clean
