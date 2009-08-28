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
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/imc.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/csg_stat.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lgmx -lboost_program_options -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../../src/tools/csg_stat

../../src/tools/csg_stat: ${BUILD_SUBPROJECTS}

../../src/tools/csg_stat: ${OBJECTFILES}
	${MKDIR} -p ../../src/tools
	${LINK.cc} -o ../../src/tools/csg_stat ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/imc.o: ../../src/tools/imc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/imc.o ../../src/tools/imc.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/csg_stat.o: ../../src/tools/csg_stat.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/src/csg/src/tools/csg_stat.o ../../src/tools/csg_stat.cc

# Subprojects
.build-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../src/tools/csg_stat

# Subprojects
.clean-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug clean
