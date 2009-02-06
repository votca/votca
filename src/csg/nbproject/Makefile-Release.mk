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
OBJECTDIR=build/Release/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/tabulatedpotential.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/bondedstatistics.o \
	${OBJECTDIR}/stdanalysis.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} dist/Release/GNU-Linux-x86/csg

dist/Release/GNU-Linux-x86/csg: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o dist/Release/GNU-Linux-x86/csg ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/tabulatedpotential.o: tabulatedpotential.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/tabulatedpotential.o tabulatedpotential.cc

${OBJECTDIR}/main.o: main.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/main.o main.cc

${OBJECTDIR}/bondedstatistics.o: bondedstatistics.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/bondedstatistics.o bondedstatistics.cc

${OBJECTDIR}/stdanalysis.o: stdanalysis.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/stdanalysis.o stdanalysis.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/csg

# Subprojects
.clean-subprojects:
