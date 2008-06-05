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
FC=

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/tmp/svntest/csg/netbeans/csg_nemat/../../src/csg_nemat.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../../../libcsg/dist/libcsg.a -lgmx -lboost_program_options -lfftw3 -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../bin/csg_nemat

../../bin/csg_nemat: ${BUILD_SUBPROJECTS}

../../bin/csg_nemat: ${OBJECTFILES}
	${MKDIR} -p ../../bin
	${LINK.cc} -o ../../bin/csg_nemat ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/tmp/svntest/csg/netbeans/csg_nemat/../../src/csg_nemat.o: ../../src/csg_nemat.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/tmp/svntest/csg/netbeans/csg_nemat/../../src
	$(COMPILE.cc) -g -I../../../libcsg -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/tmp/svntest/csg/netbeans/csg_nemat/../../src/csg_nemat.o ../../src/csg_nemat.cc

# Subprojects
.build-subprojects:
	cd ../../../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../bin/csg_nemat

# Subprojects
.clean-subprojects:
	cd ../../../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
