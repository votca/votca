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
CND_CONF=Debug
CND_DISTDIR=dist

# Include project Makefile
include Makefile_nb

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/bondedstatistics.o \
	${OBJECTDIR}/stdanalysis.o \
	${OBJECTDIR}/tabulatedpotential.o


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
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../../netbeans/libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lgmx -lexpat -lm -lfftw3 -lboost_program_options

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Debug.mk csg_boltzmann

csg_boltzmann: ../../netbeans/libcsg/../../src/libcsg/libcsg.a

csg_boltzmann: ../../../tools/netbeans/libtools/../../src/libtools/libtools.a

csg_boltzmann: ${OBJECTFILES}
	${LINK.cc} -o csg_boltzmann ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: main.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cc

${OBJECTDIR}/bondedstatistics.o: bondedstatistics.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bondedstatistics.o bondedstatistics.cc

${OBJECTDIR}/stdanalysis.o: stdanalysis.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/stdanalysis.o stdanalysis.cc

${OBJECTDIR}/tabulatedpotential.o: tabulatedpotential.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/tabulatedpotential.o tabulatedpotential.cc

# Subprojects
.build-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} csg_boltzmann

# Subprojects
.clean-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
