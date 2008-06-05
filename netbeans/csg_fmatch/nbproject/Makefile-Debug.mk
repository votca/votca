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
	${OBJECTDIR}/csg_fmatch.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../../../libcsg/dist/libcsg.a -lboost_program_options -lgmx -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../bin/csg_fmatch

../../bin/csg_fmatch: ${BUILD_SUBPROJECTS}

../../bin/csg_fmatch: ${OBJECTFILES}
	${MKDIR} -p ../../bin
	${LINK.cc} -o ../../bin/csg_fmatch ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/csg_fmatch.o: csg_fmatch.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -I../../../libcsg -o ${OBJECTDIR}/csg_fmatch.o csg_fmatch.cc

# Subprojects
.build-subprojects:
	cd ../../../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../bin/csg_fmatch

# Subprojects
.clean-subprojects:
	cd ../../../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
