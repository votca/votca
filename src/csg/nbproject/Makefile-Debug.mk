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
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../../netbeans/libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lgmx -lxml2 -lm -lfftw3 -lboost_program_options

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} ../bin/csg

../bin/csg: ${BUILD_SUBPROJECTS}

../bin/csg: ${OBJECTFILES}
	${MKDIR} -p ../bin
	${LINK.cc} -o ../bin/csg ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/tabulatedpotential.o: tabulatedpotential.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -o ${OBJECTDIR}/tabulatedpotential.o tabulatedpotential.cc

${OBJECTDIR}/main.o: main.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -o ${OBJECTDIR}/main.o main.cc

${OBJECTDIR}/bondedstatistics.o: bondedstatistics.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -o ${OBJECTDIR}/bondedstatistics.o bondedstatistics.cc

${OBJECTDIR}/stdanalysis.o: stdanalysis.cc 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -o ${OBJECTDIR}/stdanalysis.o stdanalysis.cc

# Subprojects
.build-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../bin/csg

# Subprojects
.clean-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug clean
