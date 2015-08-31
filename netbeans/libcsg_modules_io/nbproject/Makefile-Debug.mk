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
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1332856960/dlpolytopologyreader.o \
	${OBJECTDIR}/_ext/1332856960/dlpolytrajectoryreader.o \
	${OBJECTDIR}/_ext/1332856960/dlpolytrajectorywriter.o \
	${OBJECTDIR}/_ext/1332856960/gmx_print_version.o \
	${OBJECTDIR}/_ext/1332856960/gmx_version.o \
	${OBJECTDIR}/_ext/1332856960/gmx_version_check.o \
	${OBJECTDIR}/_ext/1332856960/gmx_version_nb.o \
	${OBJECTDIR}/_ext/1332856960/gmxtopologyreader.o \
	${OBJECTDIR}/_ext/1332856960/gmxtrajectoryreader.o \
	${OBJECTDIR}/_ext/1332856960/gmxtrajectorywriter.o \
	${OBJECTDIR}/_ext/1332856960/groreader.o \
	${OBJECTDIR}/_ext/1332856960/growriter.o \
	${OBJECTDIR}/_ext/1332856960/h5mdtrajectoryreader.o \
	${OBJECTDIR}/_ext/1332856960/lammpsreader.o \
	${OBJECTDIR}/_ext/1332856960/mdptopologyreader.o \
	${OBJECTDIR}/_ext/1332856960/mdptrajectoryreader.o \
	${OBJECTDIR}/_ext/1332856960/pdbreader.o \
	${OBJECTDIR}/_ext/1332856960/pdbwriter.o \
	${OBJECTDIR}/_ext/1332856960/xmltopologyreader.o \
	${OBJECTDIR}/_ext/1332856960/xyzreader.o \
	${OBJECTDIR}/_ext/1332856960/xyzwriter.o


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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a

${OBJECTDIR}/_ext/1332856960/dlpolytopologyreader.o: ../../src/libcsg/modules/io/dlpolytopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/dlpolytopologyreader.o ../../src/libcsg/modules/io/dlpolytopologyreader.cc

${OBJECTDIR}/_ext/1332856960/dlpolytrajectoryreader.o: ../../src/libcsg/modules/io/dlpolytrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/dlpolytrajectoryreader.o ../../src/libcsg/modules/io/dlpolytrajectoryreader.cc

${OBJECTDIR}/_ext/1332856960/dlpolytrajectorywriter.o: ../../src/libcsg/modules/io/dlpolytrajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/dlpolytrajectorywriter.o ../../src/libcsg/modules/io/dlpolytrajectorywriter.cc

${OBJECTDIR}/_ext/1332856960/gmx_print_version.o: ../../src/libcsg/modules/io/gmx_print_version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmx_print_version.o ../../src/libcsg/modules/io/gmx_print_version.cc

${OBJECTDIR}/_ext/1332856960/gmx_version.o: ../../src/libcsg/modules/io/gmx_version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmx_version.o ../../src/libcsg/modules/io/gmx_version.cc

${OBJECTDIR}/_ext/1332856960/gmx_version_check.o: ../../src/libcsg/modules/io/gmx_version_check.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmx_version_check.o ../../src/libcsg/modules/io/gmx_version_check.cc

${OBJECTDIR}/_ext/1332856960/gmx_version_nb.o: ../../src/libcsg/modules/io/gmx_version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmx_version_nb.o ../../src/libcsg/modules/io/gmx_version_nb.cc

${OBJECTDIR}/_ext/1332856960/gmxtopologyreader.o: ../../src/libcsg/modules/io/gmxtopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmxtopologyreader.o ../../src/libcsg/modules/io/gmxtopologyreader.cc

${OBJECTDIR}/_ext/1332856960/gmxtrajectoryreader.o: ../../src/libcsg/modules/io/gmxtrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmxtrajectoryreader.o ../../src/libcsg/modules/io/gmxtrajectoryreader.cc

${OBJECTDIR}/_ext/1332856960/gmxtrajectorywriter.o: ../../src/libcsg/modules/io/gmxtrajectorywriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/gmxtrajectorywriter.o ../../src/libcsg/modules/io/gmxtrajectorywriter.cc

${OBJECTDIR}/_ext/1332856960/groreader.o: ../../src/libcsg/modules/io/groreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/groreader.o ../../src/libcsg/modules/io/groreader.cc

${OBJECTDIR}/_ext/1332856960/growriter.o: ../../src/libcsg/modules/io/growriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/growriter.o ../../src/libcsg/modules/io/growriter.cc

${OBJECTDIR}/_ext/1332856960/h5mdtrajectoryreader.o: ../../src/libcsg/modules/io/h5mdtrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/h5mdtrajectoryreader.o ../../src/libcsg/modules/io/h5mdtrajectoryreader.cc

${OBJECTDIR}/_ext/1332856960/lammpsreader.o: ../../src/libcsg/modules/io/lammpsreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/lammpsreader.o ../../src/libcsg/modules/io/lammpsreader.cc

${OBJECTDIR}/_ext/1332856960/mdptopologyreader.o: ../../src/libcsg/modules/io/mdptopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/mdptopologyreader.o ../../src/libcsg/modules/io/mdptopologyreader.cc

${OBJECTDIR}/_ext/1332856960/mdptrajectoryreader.o: ../../src/libcsg/modules/io/mdptrajectoryreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/mdptrajectoryreader.o ../../src/libcsg/modules/io/mdptrajectoryreader.cc

${OBJECTDIR}/_ext/1332856960/pdbreader.o: ../../src/libcsg/modules/io/pdbreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/pdbreader.o ../../src/libcsg/modules/io/pdbreader.cc

${OBJECTDIR}/_ext/1332856960/pdbwriter.o: ../../src/libcsg/modules/io/pdbwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/pdbwriter.o ../../src/libcsg/modules/io/pdbwriter.cc

${OBJECTDIR}/_ext/1332856960/xmltopologyreader.o: ../../src/libcsg/modules/io/xmltopologyreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/xmltopologyreader.o ../../src/libcsg/modules/io/xmltopologyreader.cc

${OBJECTDIR}/_ext/1332856960/xyzreader.o: ../../src/libcsg/modules/io/xyzreader.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/xyzreader.o ../../src/libcsg/modules/io/xyzreader.cc

${OBJECTDIR}/_ext/1332856960/xyzwriter.o: ../../src/libcsg/modules/io/xyzwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1332856960
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1332856960/xyzwriter.o ../../src/libcsg/modules/io/xyzwriter.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibcsg_modules_io.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
