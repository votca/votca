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
FC=

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/Release/GNU-Linux-x86

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../stdanalysis.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../tabulatedpotential.o \
	${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../main.o

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
.build-conf: ${BUILD_SUBPROJECTS} dist/Release/GNU-Linux-x86/netbeans

dist/Release/GNU-Linux-x86/netbeans: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o dist/Release/GNU-Linux-x86/netbeans ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../stdanalysis.o: ../stdanalysis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/..
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../stdanalysis.o ../stdanalysis.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../tabulatedpotential.o: ../tabulatedpotential.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/..
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../tabulatedpotential.o ../tabulatedpotential.cc

${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../main.o: ../main.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/..
	$(COMPILE.cc) -O2 -o ${OBJECTDIR}/_ext/people/thnfs/homes/ruehle/svn/csg/csg/netbeans/../main.o ../main.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/netbeans

# Subprojects
.clean-subprojects:
