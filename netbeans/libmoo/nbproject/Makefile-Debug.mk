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
	${OBJECTDIR}/_ext/484467347/orbitals.o \
	${OBJECTDIR}/_ext/484467347/crgunittype.o \
	${OBJECTDIR}/_ext/484467347/crgunit.o \
	${OBJECTDIR}/_ext/484467347/charges.o \
	${OBJECTDIR}/_ext/484467347/fock.o \
	${OBJECTDIR}/_ext/484467347/basis_set.o \
	${OBJECTDIR}/_ext/484467347/jcalc2.o \
	${OBJECTDIR}/_ext/484467347/mol_and_orb.o \
	${OBJECTDIR}/_ext/484467347/units.o \
	${OBJECTDIR}/_ext/484467347/jcalc.o


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
	"${MAKE}"  -f nbproject/Makefile-Debug.mk ../../src/libmoo/libmoo.a

../../src/libmoo/libmoo.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libmoo
	${RM} ../../src/libmoo/libmoo.a
	${AR} -rv ../../src/libmoo/libmoo.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libmoo/libmoo.a

${OBJECTDIR}/_ext/484467347/orbitals.o: ../../src/libmoo/orbitals.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/orbitals.o ../../src/libmoo/orbitals.cpp

${OBJECTDIR}/_ext/484467347/crgunittype.o: ../../src/libmoo/crgunittype.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/crgunittype.o ../../src/libmoo/crgunittype.cc

${OBJECTDIR}/_ext/484467347/crgunit.o: ../../src/libmoo/crgunit.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/crgunit.o ../../src/libmoo/crgunit.cc

${OBJECTDIR}/_ext/484467347/charges.o: ../../src/libmoo/charges.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/charges.o ../../src/libmoo/charges.cpp

${OBJECTDIR}/_ext/484467347/fock.o: ../../src/libmoo/fock.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/fock.o ../../src/libmoo/fock.cc

${OBJECTDIR}/_ext/484467347/basis_set.o: ../../src/libmoo/basis_set.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/basis_set.o ../../src/libmoo/basis_set.cpp

${OBJECTDIR}/_ext/484467347/jcalc2.o: ../../src/libmoo/jcalc2.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/jcalc2.o ../../src/libmoo/jcalc2.cc

${OBJECTDIR}/_ext/484467347/mol_and_orb.o: ../../src/libmoo/mol_and_orb.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/mol_and_orb.o ../../src/libmoo/mol_and_orb.cpp

${OBJECTDIR}/_ext/484467347/units.o: ../../src/libmoo/units.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/units.o ../../src/libmoo/units.cc

${OBJECTDIR}/_ext/484467347/jcalc.o: ../../src/libmoo/jcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484467347
	${RM} $@.d
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484467347/jcalc.o ../../src/libmoo/jcalc.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} ../../src/libmoo/libmoo.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
