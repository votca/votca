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
include Makefile

# Object Directory
OBJECTDIR=build/Release/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg_noconf/netbeans/csg_nemat/../../src/tools/csg_nemat.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../../../../gmx/install/lib ../../../libcsg/trunk/netbeans/../dist/libcsg.a -lgmx -lboost_program_options -lfftw3 -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} dist/Release/GNU-Linux-x86/csg_nemat

dist/Release/GNU-Linux-x86/csg_nemat: ${BUILD_SUBPROJECTS}

dist/Release/GNU-Linux-x86/csg_nemat: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o dist/Release/GNU-Linux-x86/csg_nemat ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg_noconf/netbeans/csg_nemat/../../src/tools/csg_nemat.o: ../../src/tools/csg_nemat.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg_noconf/netbeans/csg_nemat/../../src/tools
	$(COMPILE.cc) -O2 -I../../../libcsg/trunk -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg_noconf/netbeans/csg_nemat/../../src/tools/csg_nemat.o ../../src/tools/csg_nemat.cc

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
