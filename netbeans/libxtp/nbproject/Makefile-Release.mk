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
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1ce03da5/ERIs.o \
	${OBJECTDIR}/_ext/1ce03da5/aobasis.o \
	${OBJECTDIR}/_ext/a30a7232/aocoulomb_g.o \
	${OBJECTDIR}/_ext/a30a7232/aodipole.o \
	${OBJECTDIR}/_ext/a30a7232/aoecp.o \
	${OBJECTDIR}/_ext/a30a7232/aoesp.o \
	${OBJECTDIR}/_ext/a30a7232/aokinetic.o \
	${OBJECTDIR}/_ext/a30a7232/aomatrix.o \
	${OBJECTDIR}/_ext/a30a7232/aomomentum.o \
	${OBJECTDIR}/_ext/a30a7232/aooverlap.o \
	${OBJECTDIR}/_ext/1ce03da5/aoshell.o \
	${OBJECTDIR}/_ext/1ce03da5/apolarsite.o \
	${OBJECTDIR}/_ext/1ce03da5/bsecoupling.o \
	${OBJECTDIR}/_ext/1ce03da5/calculatorfactory.o \
	${OBJECTDIR}/_ext/f1c74047/eoutersphere.o \
	${OBJECTDIR}/_ext/f1c74047/jobwriter.o \
	${OBJECTDIR}/_ext/1ce03da5/xtpapplication.o \
	${OBJECTDIR}/_ext/cdb9090a/dftengine.o \
	${OBJECTDIR}/_ext/1ce03da5/espfit.o \
	${OBJECTDIR}/_ext/1ce03da5/ewaldactor.o \
	${OBJECTDIR}/_ext/1ce03da5/extractorfactory.o \
	${OBJECTDIR}/_ext/1ce03da5/fragment.o \
	${OBJECTDIR}/_ext/1ce03da5/gdma.o \
	${OBJECTDIR}/_ext/1ce03da5/grid.o \
	${OBJECTDIR}/_ext/1ce03da5/gsl_boost_ublas_matrix_prod.o \
	${OBJECTDIR}/_ext/3f797e5a/bse.o \
	${OBJECTDIR}/_ext/3f797e5a/gwa.o \
	${OBJECTDIR}/_ext/3f797e5a/gwbse.o \
	${OBJECTDIR}/_ext/3f797e5a/rpa.o \
	${OBJECTDIR}/_ext/1ce03da5/job.o \
	${OBJECTDIR}/_ext/1ce03da5/jobapplication.o \
	${OBJECTDIR}/_ext/1ce03da5/jobcalculatorfactory.o \
	${OBJECTDIR}/_ext/d63b377e/egwbse.o \
	${OBJECTDIR}/_ext/d63b377e/idft.o \
	${OBJECTDIR}/_ext/d63b377e/iexcitoncl.o \
	${OBJECTDIR}/_ext/d63b377e/igwbse.o \
	${OBJECTDIR}/_ext/1ce03da5/molecule.o \
	${OBJECTDIR}/_ext/752fbf23/numerical_integrations.o \
	${OBJECTDIR}/_ext/752fbf23/radial_euler_maclaurin_rule.o \
	${OBJECTDIR}/_ext/752fbf23/sphere_lebedev_rule.o \
	${OBJECTDIR}/_ext/1ce03da5/orbitals.o \
	${OBJECTDIR}/_ext/1ce03da5/overlap.o \
	${OBJECTDIR}/_ext/1ce03da5/parallelpaircalc.o \
	${OBJECTDIR}/_ext/1ce03da5/parallelxjobcalc.o \
	${OBJECTDIR}/_ext/1ce03da5/poissongrid.o \
	${OBJECTDIR}/_ext/1ce03da5/polarbackground.o \
	${OBJECTDIR}/_ext/1ce03da5/polarfrag.o \
	${OBJECTDIR}/_ext/1ce03da5/polarseg.o \
	${OBJECTDIR}/_ext/1ce03da5/polarsite.o \
	${OBJECTDIR}/_ext/1ce03da5/polartop.o \
	${OBJECTDIR}/_ext/1ce03da5/progressobserver.o \
	${OBJECTDIR}/_ext/1ce03da5/qmapemachine.o \
	${OBJECTDIR}/_ext/1ce03da5/qmcalculator.o \
	${OBJECTDIR}/_ext/1ce03da5/qmdatabase.o \
	${OBJECTDIR}/_ext/1ce03da5/qmmachine.o \
	${OBJECTDIR}/_ext/1ce03da5/qmnblist.o \
	${OBJECTDIR}/_ext/1ce03da5/qmpackagefactory.o \
	${OBJECTDIR}/_ext/d95391b3/gaussian.o \
	${OBJECTDIR}/_ext/d95391b3/nwchem.o \
	${OBJECTDIR}/_ext/d95391b3/orca.o \
	${OBJECTDIR}/_ext/d95391b3/turbomole.o \
	${OBJECTDIR}/_ext/1ce03da5/qmpair.o \
	${OBJECTDIR}/_ext/1ce03da5/qmtool.o \
	${OBJECTDIR}/_ext/1ce03da5/segment.o \
	${OBJECTDIR}/_ext/1ce03da5/segmenttype.o \
	${OBJECTDIR}/_ext/1ce03da5/sqlapplication.o \
	${OBJECTDIR}/_ext/1ce03da5/statesaversqlite.o \
	${OBJECTDIR}/_ext/1ce03da5/threecenter_rep.o \
	${OBJECTDIR}/_ext/1ce03da5/threecenters.o \
	${OBJECTDIR}/_ext/1ce03da5/threecenters_dft.o \
	${OBJECTDIR}/_ext/1ce03da5/threecenters_tools.o \
	${OBJECTDIR}/_ext/1ce03da5/toolfactory.o \
	${OBJECTDIR}/_ext/402d3cf1/molpol.o \
	${OBJECTDIR}/_ext/1ce03da5/topology.o \
	${OBJECTDIR}/_ext/1ce03da5/version.o \
	${OBJECTDIR}/_ext/1ce03da5/version_nb.o \
	${OBJECTDIR}/_ext/1ce03da5/xinductor.o \
	${OBJECTDIR}/_ext/1ce03da5/xinteractor.o \
	${OBJECTDIR}/_ext/1ce03da5/xjob.o \
	${OBJECTDIR}/_ext/1ce03da5/xmapper.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a

${OBJECTDIR}/_ext/1ce03da5/ERIs.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/ERIs.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/ERIs.o ../../src/libxtp/ERIs.cc

${OBJECTDIR}/_ext/1ce03da5/aobasis.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aobasis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/aobasis.o ../../src/libxtp/aobasis.cc

${OBJECTDIR}/_ext/a30a7232/aocoulomb_g.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aocoulomb_g.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aocoulomb_g.o ../../src/libxtp/aomatrices/aocoulomb_g.cc

${OBJECTDIR}/_ext/a30a7232/aodipole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aodipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aodipole.o ../../src/libxtp/aomatrices/aodipole.cc

${OBJECTDIR}/_ext/a30a7232/aoecp.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aoecp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aoecp.o ../../src/libxtp/aomatrices/aoecp.cc

${OBJECTDIR}/_ext/a30a7232/aoesp.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aoesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aoesp.o ../../src/libxtp/aomatrices/aoesp.cc

${OBJECTDIR}/_ext/a30a7232/aokinetic.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aokinetic.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aokinetic.o ../../src/libxtp/aomatrices/aokinetic.cc

${OBJECTDIR}/_ext/a30a7232/aomatrix.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aomatrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aomatrix.o ../../src/libxtp/aomatrices/aomatrix.cc

${OBJECTDIR}/_ext/a30a7232/aomomentum.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aomomentum.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aomomentum.o ../../src/libxtp/aomatrices/aomomentum.cc

${OBJECTDIR}/_ext/a30a7232/aooverlap.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aooverlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/a30a7232
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/a30a7232/aooverlap.o ../../src/libxtp/aomatrices/aooverlap.cc

${OBJECTDIR}/_ext/1ce03da5/aoshell.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aoshell.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/aoshell.o ../../src/libxtp/aoshell.cc

${OBJECTDIR}/_ext/1ce03da5/apolarsite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/apolarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/apolarsite.o ../../src/libxtp/apolarsite.cc

${OBJECTDIR}/_ext/1ce03da5/bsecoupling.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/bsecoupling.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/bsecoupling.o ../../src/libxtp/bsecoupling.cc

${OBJECTDIR}/_ext/1ce03da5/calculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/calculatorfactory.o ../../src/libxtp/calculatorfactory.cc

${OBJECTDIR}/_ext/f1c74047/eoutersphere.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f1c74047
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f1c74047/eoutersphere.o ../../src/libxtp/calculators/eoutersphere.cc

${OBJECTDIR}/_ext/f1c74047/jobwriter.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/calculators/jobwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f1c74047
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f1c74047/jobwriter.o ../../src/libxtp/calculators/jobwriter.cc

${OBJECTDIR}/_ext/1ce03da5/xtpapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xtpapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/xtpapplication.o ../../src/libxtp/xtpapplication.cc

${OBJECTDIR}/_ext/cdb9090a/dftengine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/dftengine/dftengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/cdb9090a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/cdb9090a/dftengine.o ../../src/libxtp/dftengine/dftengine.cc

${OBJECTDIR}/_ext/1ce03da5/espfit.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/espfit.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/espfit.o ../../src/libxtp/espfit.cc

${OBJECTDIR}/_ext/1ce03da5/ewaldactor.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/ewaldactor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/ewaldactor.o ../../src/libxtp/ewaldactor.cc

${OBJECTDIR}/_ext/1ce03da5/extractorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/extractorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/extractorfactory.o ../../src/libxtp/extractorfactory.cc

${OBJECTDIR}/_ext/1ce03da5/fragment.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/fragment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/fragment.o ../../src/libxtp/fragment.cc

${OBJECTDIR}/_ext/1ce03da5/gdma.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gdma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/gdma.o ../../src/libxtp/gdma.cc

${OBJECTDIR}/_ext/1ce03da5/grid.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/grid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/grid.o ../../src/libxtp/grid.cc

${OBJECTDIR}/_ext/1ce03da5/gsl_boost_ublas_matrix_prod.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gsl_boost_ublas_matrix_prod.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/gsl_boost_ublas_matrix_prod.o ../../src/libxtp/gsl_boost_ublas_matrix_prod.cc

${OBJECTDIR}/_ext/3f797e5a/bse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/bse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/3f797e5a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/3f797e5a/bse.o ../../src/libxtp/gwbse/bse.cc

${OBJECTDIR}/_ext/3f797e5a/gwa.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/gwa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/3f797e5a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/3f797e5a/gwa.o ../../src/libxtp/gwbse/gwa.cc

${OBJECTDIR}/_ext/3f797e5a/gwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/gwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/3f797e5a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/3f797e5a/gwbse.o ../../src/libxtp/gwbse/gwbse.cc

${OBJECTDIR}/_ext/3f797e5a/rpa.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/rpa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/3f797e5a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/3f797e5a/rpa.o ../../src/libxtp/gwbse/rpa.cc

${OBJECTDIR}/_ext/1ce03da5/job.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/job.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/job.o ../../src/libxtp/job.cc

${OBJECTDIR}/_ext/1ce03da5/jobapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/jobapplication.o ../../src/libxtp/jobapplication.cc

${OBJECTDIR}/_ext/1ce03da5/jobcalculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/jobcalculatorfactory.o ../../src/libxtp/jobcalculatorfactory.cc

${OBJECTDIR}/_ext/d63b377e/egwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/egwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d63b377e
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d63b377e/egwbse.o ../../src/libxtp/jobcalculators/egwbse.cc

${OBJECTDIR}/_ext/d63b377e/idft.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d63b377e
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d63b377e/idft.o ../../src/libxtp/jobcalculators/idft.cc

${OBJECTDIR}/_ext/d63b377e/iexcitoncl.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/iexcitoncl.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d63b377e
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d63b377e/iexcitoncl.o ../../src/libxtp/jobcalculators/iexcitoncl.cc

${OBJECTDIR}/_ext/d63b377e/igwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/igwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d63b377e
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d63b377e/igwbse.o ../../src/libxtp/jobcalculators/igwbse.cc

${OBJECTDIR}/_ext/1ce03da5/molecule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/molecule.o ../../src/libxtp/molecule.cc

${OBJECTDIR}/_ext/752fbf23/numerical_integrations.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/numerical_integrations.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/752fbf23
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/752fbf23/numerical_integrations.o ../../src/libxtp/numerical_integration/numerical_integrations.cc

${OBJECTDIR}/_ext/752fbf23/radial_euler_maclaurin_rule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/752fbf23
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/752fbf23/radial_euler_maclaurin_rule.o ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc

${OBJECTDIR}/_ext/752fbf23/sphere_lebedev_rule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/752fbf23
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/752fbf23/sphere_lebedev_rule.o ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc

${OBJECTDIR}/_ext/1ce03da5/orbitals.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/orbitals.o ../../src/libxtp/orbitals.cc

${OBJECTDIR}/_ext/1ce03da5/overlap.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/overlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/overlap.o ../../src/libxtp/overlap.cc

${OBJECTDIR}/_ext/1ce03da5/parallelpaircalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/parallelpaircalc.o ../../src/libxtp/parallelpaircalc.cc

${OBJECTDIR}/_ext/1ce03da5/parallelxjobcalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/parallelxjobcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/parallelxjobcalc.o ../../src/libxtp/parallelxjobcalc.cc

${OBJECTDIR}/_ext/1ce03da5/poissongrid.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/poissongrid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/poissongrid.o ../../src/libxtp/poissongrid.cc

${OBJECTDIR}/_ext/1ce03da5/polarbackground.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/polarbackground.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/polarbackground.o ../../src/libxtp/polarbackground.cc

${OBJECTDIR}/_ext/1ce03da5/polarfrag.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/polarfrag.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/polarfrag.o ../../src/libxtp/polarfrag.cc

${OBJECTDIR}/_ext/1ce03da5/polarseg.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/polarseg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/polarseg.o ../../src/libxtp/polarseg.cc

${OBJECTDIR}/_ext/1ce03da5/polarsite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/polarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/polarsite.o ../../src/libxtp/polarsite.cc

${OBJECTDIR}/_ext/1ce03da5/polartop.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/polartop.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/polartop.o ../../src/libxtp/polartop.cc

${OBJECTDIR}/_ext/1ce03da5/progressobserver.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/progressobserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/progressobserver.o ../../src/libxtp/progressobserver.cc

${OBJECTDIR}/_ext/1ce03da5/qmapemachine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmapemachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmapemachine.o ../../src/libxtp/qmapemachine.cc

${OBJECTDIR}/_ext/1ce03da5/qmcalculator.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmcalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmcalculator.o ../../src/libxtp/qmcalculator.cc

${OBJECTDIR}/_ext/1ce03da5/qmdatabase.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmdatabase.o ../../src/libxtp/qmdatabase.cc

${OBJECTDIR}/_ext/1ce03da5/qmmachine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmmachine.o ../../src/libxtp/qmmachine.cc

${OBJECTDIR}/_ext/1ce03da5/qmnblist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmnblist.o ../../src/libxtp/qmnblist.cc

${OBJECTDIR}/_ext/1ce03da5/qmpackagefactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmpackagefactory.o ../../src/libxtp/qmpackagefactory.cc

${OBJECTDIR}/_ext/d95391b3/gaussian.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d95391b3
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d95391b3/gaussian.o ../../src/libxtp/qmpackages/gaussian.cc

${OBJECTDIR}/_ext/d95391b3/nwchem.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/nwchem.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d95391b3
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d95391b3/nwchem.o ../../src/libxtp/qmpackages/nwchem.cc

${OBJECTDIR}/_ext/d95391b3/orca.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/orca.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d95391b3
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d95391b3/orca.o ../../src/libxtp/qmpackages/orca.cc

${OBJECTDIR}/_ext/d95391b3/turbomole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/turbomole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/d95391b3
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/d95391b3/turbomole.o ../../src/libxtp/qmpackages/turbomole.cc

${OBJECTDIR}/_ext/1ce03da5/qmpair.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmpair.o ../../src/libxtp/qmpair.cc

${OBJECTDIR}/_ext/1ce03da5/qmtool.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmtool.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/qmtool.o ../../src/libxtp/qmtool.cc

${OBJECTDIR}/_ext/1ce03da5/segment.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/segment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/segment.o ../../src/libxtp/segment.cc

${OBJECTDIR}/_ext/1ce03da5/segmenttype.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/segmenttype.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/segmenttype.o ../../src/libxtp/segmenttype.cc

${OBJECTDIR}/_ext/1ce03da5/sqlapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/sqlapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/sqlapplication.o ../../src/libxtp/sqlapplication.cc

${OBJECTDIR}/_ext/1ce03da5/statesaversqlite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/statesaversqlite.o ../../src/libxtp/statesaversqlite.cc

${OBJECTDIR}/_ext/1ce03da5/threecenter_rep.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenter_rep.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/threecenter_rep.o ../../src/libxtp/threecenter_rep.cc

${OBJECTDIR}/_ext/1ce03da5/threecenters.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/threecenters.o ../../src/libxtp/threecenters.cc

${OBJECTDIR}/_ext/1ce03da5/threecenters_dft.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters_dft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/threecenters_dft.o ../../src/libxtp/threecenters_dft.cc

${OBJECTDIR}/_ext/1ce03da5/threecenters_tools.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters_tools.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/threecenters_tools.o ../../src/libxtp/threecenters_tools.cc

${OBJECTDIR}/_ext/1ce03da5/toolfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/toolfactory.o ../../src/libxtp/toolfactory.cc

${OBJECTDIR}/_ext/402d3cf1/molpol.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/tools/molpol.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/402d3cf1
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/402d3cf1/molpol.o ../../src/libxtp/tools/molpol.cc

${OBJECTDIR}/_ext/1ce03da5/topology.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/topology.o ../../src/libxtp/topology.cc

${OBJECTDIR}/_ext/1ce03da5/version.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/version.o ../../src/libxtp/version.cc

${OBJECTDIR}/_ext/1ce03da5/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/version_nb.o ../../src/libxtp/version_nb.cc

${OBJECTDIR}/_ext/1ce03da5/xinductor.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xinductor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/xinductor.o ../../src/libxtp/xinductor.cc

${OBJECTDIR}/_ext/1ce03da5/xinteractor.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xinteractor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/xinteractor.o ../../src/libxtp/xinteractor.cc

${OBJECTDIR}/_ext/1ce03da5/xjob.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xjob.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/xjob.o ../../src/libxtp/xjob.cc

${OBJECTDIR}/_ext/1ce03da5/xmapper.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xmapper.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce03da5
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce03da5/xmapper.o ../../src/libxtp/xmapper.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibxtp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
