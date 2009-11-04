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

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/people/thnfs/homes/ruehle/gmx/lib ../../netbeans/libcsg/../../src/libcsg/libcsg.a ../../../tools/netbeans/libtools/../../src/libtools/libtools.a -lgmx -lxml2 -lm -lfftw3 -lboost_program_options

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk csg

csg: ../../netbeans/libcsg/../../src/libcsg/libcsg.a

csg: ../../../tools/netbeans/libtools/../../src/libtools/libtools.a

csg: ${OBJECTFILES}
	${LINK.cc} -o csg ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/tabulatedpotential.o: nbproject/Makefile-${CND_CONF}.mk tabulatedpotential.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/tabulatedpotential.o tabulatedpotential.cc

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cc

${OBJECTDIR}/bondedstatistics.o: nbproject/Makefile-${CND_CONF}.mk bondedstatistics.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/bondedstatistics.o bondedstatistics.cc

${OBJECTDIR}/stdanalysis.o: nbproject/Makefile-${CND_CONF}.mk stdanalysis.cc 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/people/thnfs/homes/ruehle/gmx/install/include/gromacs -I../../include -I../../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/stdanalysis.o stdanalysis.cc

# Subprojects
.build-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} csg

# Subprojects
.clean-subprojects:
	cd ../../netbeans/libcsg && ${MAKE}  -f Makefile_nb CONF=Debug clean
	cd ../../../tools/netbeans/libtools && ${MAKE}  -f Makefile_nb CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
