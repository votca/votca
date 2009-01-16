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
FC=g95

# Include project Makefile
include Makefile_nb

# Object Directory
OBJECTDIR=build/Debug/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/qm_molecule.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/basis_set.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/orbitals.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/charges.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/fock_matrix.o

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
.build-conf: ${BUILD_SUBPROJECTS} ../../src/libmoo/libmoo.a

../../src/libmoo/libmoo.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libmoo
	${RM} ../../src/libmoo/libmoo.a
	${AR} rv ../../src/libmoo/libmoo.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libmoo/libmoo.a

${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/qm_molecule.o: ../../src/libmoo/qm_molecule.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/qm_molecule.o ../../src/libmoo/qm_molecule.cpp

${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/basis_set.o: ../../src/libmoo/basis_set.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/basis_set.o ../../src/libmoo/basis_set.cpp

${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/orbitals.o: ../../src/libmoo/orbitals.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/orbitals.o ../../src/libmoo/orbitals.cpp

${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/charges.o: ../../src/libmoo/charges.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/charges.o ../../src/libmoo/charges.cpp

${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/fock_matrix.o: ../../src/libmoo/fock_matrix.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo
	$(COMPILE.cc) -g -I../../include -I../../../include -o ${OBJECTDIR}/_ext/people/thnfs/homes/vehoff/src/csg_svn/moo/src/libmoo/fock_matrix.o ../../src/libmoo/fock_matrix.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} ../../src/libmoo/libmoo.a

# Subprojects
.clean-subprojects:
