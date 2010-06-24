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
FC=g95
AS=

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
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/cubicspline.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/tokenizer.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parcer.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/property.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/random.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogramnew.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/matrix.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/crosscorrelate.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/correlate.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/rangeparser.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogram.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parsexml.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/table.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/linalg.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/version_nb.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/datacollection.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-I../../include
CXXFLAGS=-I../../include

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk ../../src/libtools/libtools.a

../../src/libtools/libtools.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libtools
	${RM} ../../src/libtools/libtools.a
	${AR} rv ../../src/libtools/libtools.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libtools/libtools.a

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/cubicspline.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/cubicspline.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/cubicspline.o ../../src/libtools/cubicspline.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/tokenizer.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/tokenizer.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/tokenizer.o ../../src/libtools/tokenizer.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parcer.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/parcer.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parcer.o ../../src/libtools/parcer.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/property.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/property.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/property.o ../../src/libtools/property.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/random.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/random.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/random.o ../../src/libtools/random.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogramnew.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/histogramnew.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogramnew.o ../../src/libtools/histogramnew.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/matrix.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/matrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/matrix.o ../../src/libtools/matrix.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/crosscorrelate.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/crosscorrelate.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/crosscorrelate.o ../../src/libtools/crosscorrelate.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/correlate.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/correlate.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/correlate.o ../../src/libtools/correlate.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/rangeparser.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/rangeparser.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/rangeparser.o ../../src/libtools/rangeparser.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogram.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/histogram.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/histogram.o ../../src/libtools/histogram.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parsexml.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/parsexml.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/parsexml.o ../../src/libtools/parsexml.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/table.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/table.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/table.o ../../src/libtools/table.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/linalg.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/linalg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/linalg.o ../../src/libtools/linalg.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/version_nb.o ../../src/libtools/version_nb.cc

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/datacollection.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libtools/datacollection.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/src/libtools/datacollection.o ../../src/libtools/datacollection.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} ../../src/libtools/libtools.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
