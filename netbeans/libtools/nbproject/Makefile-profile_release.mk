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
include Makefile_nb

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1723120609/version_nb.o \
	${OBJECTDIR}/_ext/1723120609/easylock.o \
	${OBJECTDIR}/_ext/1723120609/cubicspline.o \
	${OBJECTDIR}/_ext/1723120609/application.o \
	${OBJECTDIR}/_ext/1723120609/tokenizer.o \
	${OBJECTDIR}/_ext/1723120609/table.o \
	${OBJECTDIR}/_ext/1723120609/parcer.o \
	${OBJECTDIR}/_ext/1723120609/random.o \
	${OBJECTDIR}/_ext/1723120609/mutex.o \
	${OBJECTDIR}/_ext/1723120609/property.o \
	${OBJECTDIR}/_ext/1723120609/akimaspline.o \
	${OBJECTDIR}/_ext/1723120609/parsexml.o \
	${OBJECTDIR}/_ext/1723120609/crosscorrelate.o \
	${OBJECTDIR}/_ext/1723120609/histogramnew.o \
	${OBJECTDIR}/_ext/1723120609/linalg.o \
	${OBJECTDIR}/_ext/1723120609/matrix.o \
	${OBJECTDIR}/_ext/1723120609/datacollection.o \
	${OBJECTDIR}/_ext/1723120609/correlate.o \
	${OBJECTDIR}/_ext/1723120609/histogram.o \
	${OBJECTDIR}/_ext/1723120609/linspline.o \
	${OBJECTDIR}/_ext/1723120609/thread.o \
	${OBJECTDIR}/_ext/1723120609/spline.o \
	${OBJECTDIR}/_ext/1723120609/rangeparser.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-I../../include -pg
CXXFLAGS=-I../../include -pg

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-profile_release.mk ../../src/libtools/libtools.a

../../src/libtools/libtools.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libtools
	${RM} ../../src/libtools/libtools.a
	${AR} -rv ../../src/libtools/libtools.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libtools/libtools.a

${OBJECTDIR}/_ext/1723120609/version_nb.o: ../../src/libtools/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/version_nb.o ../../src/libtools/version_nb.cc

${OBJECTDIR}/_ext/1723120609/easylock.o: ../../src/libtools/easylock.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/easylock.o ../../src/libtools/easylock.cc

${OBJECTDIR}/_ext/1723120609/cubicspline.o: ../../src/libtools/cubicspline.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/cubicspline.o ../../src/libtools/cubicspline.cc

${OBJECTDIR}/_ext/1723120609/application.o: ../../src/libtools/application.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/application.o ../../src/libtools/application.cc

${OBJECTDIR}/_ext/1723120609/tokenizer.o: ../../src/libtools/tokenizer.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/tokenizer.o ../../src/libtools/tokenizer.cc

${OBJECTDIR}/_ext/1723120609/table.o: ../../src/libtools/table.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/table.o ../../src/libtools/table.cc

${OBJECTDIR}/_ext/1723120609/parcer.o: ../../src/libtools/parcer.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/parcer.o ../../src/libtools/parcer.cc

${OBJECTDIR}/_ext/1723120609/random.o: ../../src/libtools/random.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/random.o ../../src/libtools/random.cc

${OBJECTDIR}/_ext/1723120609/mutex.o: ../../src/libtools/mutex.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/mutex.o ../../src/libtools/mutex.cc

${OBJECTDIR}/_ext/1723120609/property.o: ../../src/libtools/property.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/property.o ../../src/libtools/property.cc

${OBJECTDIR}/_ext/1723120609/akimaspline.o: ../../src/libtools/akimaspline.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/akimaspline.o ../../src/libtools/akimaspline.cc

${OBJECTDIR}/_ext/1723120609/parsexml.o: ../../src/libtools/parsexml.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/parsexml.o ../../src/libtools/parsexml.cc

${OBJECTDIR}/_ext/1723120609/crosscorrelate.o: ../../src/libtools/crosscorrelate.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/crosscorrelate.o ../../src/libtools/crosscorrelate.cc

${OBJECTDIR}/_ext/1723120609/histogramnew.o: ../../src/libtools/histogramnew.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/histogramnew.o ../../src/libtools/histogramnew.cc

${OBJECTDIR}/_ext/1723120609/linalg.o: ../../src/libtools/linalg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/linalg.o ../../src/libtools/linalg.cc

${OBJECTDIR}/_ext/1723120609/matrix.o: ../../src/libtools/matrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/matrix.o ../../src/libtools/matrix.cc

${OBJECTDIR}/_ext/1723120609/datacollection.o: ../../src/libtools/datacollection.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/datacollection.o ../../src/libtools/datacollection.cc

${OBJECTDIR}/_ext/1723120609/correlate.o: ../../src/libtools/correlate.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/correlate.o ../../src/libtools/correlate.cc

${OBJECTDIR}/_ext/1723120609/histogram.o: ../../src/libtools/histogram.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/histogram.o ../../src/libtools/histogram.cc

${OBJECTDIR}/_ext/1723120609/linspline.o: ../../src/libtools/linspline.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/linspline.o ../../src/libtools/linspline.cc

${OBJECTDIR}/_ext/1723120609/thread.o: ../../src/libtools/thread.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/thread.o ../../src/libtools/thread.cc

${OBJECTDIR}/_ext/1723120609/spline.o: ../../src/libtools/spline.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/spline.o ../../src/libtools/spline.cc

${OBJECTDIR}/_ext/1723120609/rangeparser.o: ../../src/libtools/rangeparser.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1723120609
	${RM} $@.d
	$(COMPILE.cc) -O3 -I../../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1723120609/rangeparser.o ../../src/libtools/rangeparser.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/profile_release
	${RM} ../../src/libtools/libtools.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
