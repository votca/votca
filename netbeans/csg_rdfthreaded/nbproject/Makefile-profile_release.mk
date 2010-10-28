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
CND_CONF=profile_release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/715944016/csg_rdfthreaded.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pg
CXXFLAGS=-pg

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lgmx -lboost_program_options -lxml2 -lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-profile_release.mk ../../src/tools/csg_stat

../../src/tools/csg_stat: ../libcsg/../../src/libcsg/libcsg.a

../../src/tools/csg_stat: ../../../tools/netbeans/libtools/../../src/libtools/libtools.a

../../src/tools/csg_stat: ${OBJECTFILES}
	${MKDIR} -p ../../src/tools
	${LINK.cc} -pg -o ../../src/tools/csg_stat ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/715944016/csg_rdfthreaded.o: ../../src/tools/csg_rdfthreaded.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/715944016
	${RM} $@.d
	$(COMPILE.cc) -g -O -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/715944016/csg_rdfthreaded.o ../../src/tools/csg_rdfthreaded.cc

# Subprojects
.build-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=profile_release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/profile_release
	${RM} ../../src/tools/csg_stat

# Subprojects
.clean-subprojects:
	cd ../libcsg && ${MAKE}  -f Makefile_nb CONF=profile_release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
