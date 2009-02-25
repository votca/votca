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
FC=g95

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/lukyanov/src/csg/src/tools/csg_fmatch.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lboost_program_options -lgmx -lxml2 -lgsl -lgslcblas -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../src/tools/csg_fmatch

../../src/tools/csg_fmatch: ${BUILD_SUBPROJECTS}

../../src/tools/csg_fmatch: ${OBJECTFILES}
	${MKDIR} -p ../../src/tools
	${LINK.cc} -o ../../src/tools/csg_fmatch ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/lukyanov/src/csg/src/tools/csg_fmatch.o: ../../src/tools/csg_fmatch.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/lukyanov/src/csg/src/tools
	$(COMPILE.cc) -g -O -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/lukyanov/src/csg/src/tools/csg_fmatch.o ../../src/tools/csg_fmatch.cc

# Subprojects
.build-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../src/tools/csg_fmatch

# Subprojects
.clean-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug clean
