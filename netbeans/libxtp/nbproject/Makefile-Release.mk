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
	${OBJECTDIR}/_ext/1ce08c7a/ERIs.o \
	${OBJECTDIR}/_ext/1ce08c7a/aobasis.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aocoulomb_g.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o \
	${OBJECTDIR}/_ext/1ce08c7a/aoshell.o \
	${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o \
	${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o \
	${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o \
	${OBJECTDIR}/_ext/5e78919f/dftengine.o \
	${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o \
	${OBJECTDIR}/_ext/1ce08c7a/espfit.o \
	${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o \
	${OBJECTDIR}/_ext/1ce08c7a/gdma.o \
	${OBJECTDIR}/_ext/1ce08c7a/grid.o \
	${OBJECTDIR}/_ext/6916996f/bse.o \
	${OBJECTDIR}/_ext/6916996f/gwa.o \
	${OBJECTDIR}/_ext/6916996f/gwbse.o \
	${OBJECTDIR}/_ext/6916996f/rpa.o \
	${OBJECTDIR}/_ext/1ce08c7a/job.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o \
	${OBJECTDIR}/_ext/f632c409/dma.o \
	${OBJECTDIR}/_ext/f632c409/egwbse.o \
	${OBJECTDIR}/_ext/f632c409/iexcitoncl.o \
	${OBJECTDIR}/_ext/f632c409/igwbse.o \
	${OBJECTDIR}/_ext/1ce08c7a/lowdin.o \
	${OBJECTDIR}/_ext/1ce08c7a/mulliken.o \
	${OBJECTDIR}/_ext/1ce08c7a/nbo.o \
	${OBJECTDIR}/_ext/4d261038/numerical_integrations.o \
	${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o \
	${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o \
	${OBJECTDIR}/_ext/1ce08c7a/orbitals.o \
	${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o \
	${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o \
	${OBJECTDIR}/_ext/60851bbe/cpmd.o \
	${OBJECTDIR}/_ext/60851bbe/gaussian.o \
	${OBJECTDIR}/_ext/60851bbe/nwchem.o \
	${OBJECTDIR}/_ext/60851bbe/orca.o \
	${OBJECTDIR}/_ext/60851bbe/turbomole.o \
	${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o \
	${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o \
	${OBJECTDIR}/_ext/1ce08c7a/version.o \
	${OBJECTDIR}/_ext/1ce08c7a/version_nb.o \
	${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o


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

${OBJECTDIR}/_ext/1ce08c7a/ERIs.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/ERIs.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/ERIs.o ../../src/libxtp/ERIs.cc

${OBJECTDIR}/_ext/1ce08c7a/aobasis.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aobasis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aobasis.o ../../src/libxtp/aobasis.cc

${OBJECTDIR}/_ext/2a3bfc3d/aocoulomb_g.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aocoulomb_g.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aocoulomb_g.o ../../src/libxtp/aomatrices/aocoulomb_g.cc

${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aodipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o ../../src/libxtp/aomatrices/aodipole.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aoecp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o ../../src/libxtp/aomatrices/aoecp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aoesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o ../../src/libxtp/aomatrices/aoesp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aokinetic.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o ../../src/libxtp/aomatrices/aokinetic.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aomatrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o ../../src/libxtp/aomatrices/aomatrix.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aomomentum.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o ../../src/libxtp/aomatrices/aomomentum.cc

${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aomatrices/aooverlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o ../../src/libxtp/aomatrices/aooverlap.cc

${OBJECTDIR}/_ext/1ce08c7a/aoshell.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/aoshell.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aoshell.o ../../src/libxtp/aoshell.cc

${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/bsecoupling.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o ../../src/libxtp/bsecoupling.cc

${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/bulkesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o ../../src/libxtp/bulkesp.cc

${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o ../../src/libxtp/calculatorfactory.cc

${OBJECTDIR}/_ext/5e78919f/dftengine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/dftengine/dftengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/5e78919f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/5e78919f/dftengine.o ../../src/libxtp/dftengine/dftengine.cc

${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/esp2multipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o ../../src/libxtp/esp2multipole.cc

${OBJECTDIR}/_ext/1ce08c7a/espfit.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/espfit.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/espfit.o ../../src/libxtp/espfit.cc

${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/extractorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o ../../src/libxtp/extractorfactory.cc

${OBJECTDIR}/_ext/1ce08c7a/gdma.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gdma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gdma.o ../../src/libxtp/gdma.cc

${OBJECTDIR}/_ext/1ce08c7a/grid.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/grid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/grid.o ../../src/libxtp/grid.cc

${OBJECTDIR}/_ext/6916996f/bse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/bse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/bse.o ../../src/libxtp/gwbse/bse.cc

${OBJECTDIR}/_ext/6916996f/gwa.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/gwa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwa.o ../../src/libxtp/gwbse/gwa.cc

${OBJECTDIR}/_ext/6916996f/gwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/gwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwbse.o ../../src/libxtp/gwbse/gwbse.cc

${OBJECTDIR}/_ext/6916996f/rpa.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/gwbse/rpa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/rpa.o ../../src/libxtp/gwbse/rpa.cc

${OBJECTDIR}/_ext/1ce08c7a/job.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/job.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/job.o ../../src/libxtp/job.cc

${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o ../../src/libxtp/jobapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o ../../src/libxtp/jobcalculatorfactory.cc

${OBJECTDIR}/_ext/f632c409/dma.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/dma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/dma.o ../../src/libxtp/jobcalculators/dma.cc

${OBJECTDIR}/_ext/f632c409/egwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/egwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/egwbse.o ../../src/libxtp/jobcalculators/egwbse.cc

${OBJECTDIR}/_ext/f632c409/iexcitoncl.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/iexcitoncl.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/iexcitoncl.o ../../src/libxtp/jobcalculators/iexcitoncl.cc

${OBJECTDIR}/_ext/f632c409/igwbse.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/jobcalculators/igwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/igwbse.o ../../src/libxtp/jobcalculators/igwbse.cc

${OBJECTDIR}/_ext/1ce08c7a/lowdin.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/lowdin.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/lowdin.o ../../src/libxtp/lowdin.cc

${OBJECTDIR}/_ext/1ce08c7a/mulliken.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/mulliken.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/mulliken.o ../../src/libxtp/mulliken.cc

${OBJECTDIR}/_ext/1ce08c7a/nbo.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/nbo.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/nbo.o ../../src/libxtp/nbo.cc

${OBJECTDIR}/_ext/4d261038/numerical_integrations.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/numerical_integrations.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/numerical_integrations.o ../../src/libxtp/numerical_integration/numerical_integrations.cc

${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc

${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc

${OBJECTDIR}/_ext/1ce08c7a/orbitals.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/orbitals.o ../../src/libxtp/orbitals.cc

${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o ../../src/libxtp/parallelpaircalc.cc

${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/parallelxjobcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o ../../src/libxtp/parallelxjobcalc.cc

${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmapemachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o ../../src/libxtp/qmapemachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmcalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o ../../src/libxtp/qmcalculator.cc

${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o ../../src/libxtp/qmdatabase.cc

${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o ../../src/libxtp/qmmachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o ../../src/libxtp/qmpackagefactory.cc

${OBJECTDIR}/_ext/60851bbe/cpmd.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/cpmd.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/cpmd.o ../../src/libxtp/qmpackages/cpmd.cc

${OBJECTDIR}/_ext/60851bbe/gaussian.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/gaussian.o ../../src/libxtp/qmpackages/gaussian.cc

${OBJECTDIR}/_ext/60851bbe/nwchem.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/nwchem.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/nwchem.o ../../src/libxtp/qmpackages/nwchem.cc

${OBJECTDIR}/_ext/60851bbe/orca.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/orca.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/orca.o ../../src/libxtp/qmpackages/orca.cc

${OBJECTDIR}/_ext/60851bbe/turbomole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/qmpackages/turbomole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/turbomole.o ../../src/libxtp/qmpackages/turbomole.cc

${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/sqlapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o ../../src/libxtp/sqlapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o ../../src/libxtp/statesaversqlite.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenter_rep.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o ../../src/libxtp/threecenter_rep.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters.o ../../src/libxtp/threecenters.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters_dft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o ../../src/libxtp/threecenters_dft.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters_g.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o ../../src/libxtp/threecenters_g.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/threecenters_tools.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o ../../src/libxtp/threecenters_tools.cc

${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o ../../src/libxtp/toolfactory.cc

${OBJECTDIR}/_ext/1ce08c7a/version.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version.o ../../src/libxtp/version.cc

${OBJECTDIR}/_ext/1ce08c7a/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version_nb.o ../../src/libxtp/version_nb.cc

${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/xtpapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o ../../src/libxtp/xtpapplication.cc

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
