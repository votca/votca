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
CND_CONF=Debug
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
	${OBJECTDIR}/_ext/1ce08c7a/apolarsite.o \
	${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o \
	${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o \
	${OBJECTDIR}/_ext/50c6f79c/eoutersphere.o \
	${OBJECTDIR}/_ext/50c6f79c/jobwriter.o \
	${OBJECTDIR}/_ext/5e78919f/dftengine.o \
	${OBJECTDIR}/_ext/1ce08c7a/espfit.o \
	${OBJECTDIR}/_ext/1ce08c7a/ewaldactor.o \
	${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o \
	${OBJECTDIR}/_ext/1ce08c7a/fragment.o \
	${OBJECTDIR}/_ext/1ce08c7a/gdma.o \
	${OBJECTDIR}/_ext/1ce08c7a/grid.o \
	${OBJECTDIR}/_ext/1ce08c7a/gsl_boost_ublas_matrix_prod.o \
	${OBJECTDIR}/_ext/6916996f/bse.o \
	${OBJECTDIR}/_ext/6916996f/gwa.o \
	${OBJECTDIR}/_ext/6916996f/gwbse.o \
	${OBJECTDIR}/_ext/6916996f/rpa.o \
	${OBJECTDIR}/_ext/1ce08c7a/job.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o \
	${OBJECTDIR}/_ext/f632c409/egwbse.o \
	${OBJECTDIR}/_ext/f632c409/idft.o \
	${OBJECTDIR}/_ext/f632c409/iexcitoncl.o \
	${OBJECTDIR}/_ext/f632c409/igwbse.o \
	${OBJECTDIR}/_ext/1ce08c7a/molecule.o \
	${OBJECTDIR}/_ext/4d261038/numerical_integrations.o \
	${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o \
	${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o \
	${OBJECTDIR}/_ext/1ce08c7a/orbitals.o \
	${OBJECTDIR}/_ext/1ce08c7a/overlap.o \
	${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o \
	${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o \
	${OBJECTDIR}/_ext/1ce08c7a/poissongrid.o \
	${OBJECTDIR}/_ext/1ce08c7a/polarbackground.o \
	${OBJECTDIR}/_ext/1ce08c7a/polarfrag.o \
	${OBJECTDIR}/_ext/1ce08c7a/polarseg.o \
	${OBJECTDIR}/_ext/1ce08c7a/polarsite.o \
	${OBJECTDIR}/_ext/1ce08c7a/polartop.o \
	${OBJECTDIR}/_ext/1ce08c7a/progressobserver.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmnblist.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o \
	${OBJECTDIR}/_ext/60851bbe/gaussian.o \
	${OBJECTDIR}/_ext/60851bbe/nwchem.o \
	${OBJECTDIR}/_ext/60851bbe/orca.o \
	${OBJECTDIR}/_ext/60851bbe/turbomole.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmpair.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmtool.o \
	${OBJECTDIR}/_ext/1ce08c7a/segment.o \
	${OBJECTDIR}/_ext/1ce08c7a/segmenttype.o \
	${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o \
	${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o \
	${OBJECTDIR}/_ext/69ca5806/molpol.o \
	${OBJECTDIR}/_ext/1ce08c7a/topology.o \
	${OBJECTDIR}/_ext/1ce08c7a/version.o \
	${OBJECTDIR}/_ext/1ce08c7a/version_nb.o \
	${OBJECTDIR}/_ext/1ce08c7a/xinductor.o \
	${OBJECTDIR}/_ext/1ce08c7a/xinteractor.o \
	${OBJECTDIR}/_ext/1ce08c7a/xjob.o \
	${OBJECTDIR}/_ext/1ce08c7a/xmapper.o \
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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/libxtp.a

../../src/libxtp/libxtp.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libxtp
	${RM} ../../src/libxtp/libxtp.a
	${AR} -rv ../../src/libxtp/libxtp.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libxtp/libxtp.a

${OBJECTDIR}/_ext/1ce08c7a/ERIs.o: ../../src/libxtp/ERIs.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/ERIs.o ../../src/libxtp/ERIs.cc

${OBJECTDIR}/_ext/1ce08c7a/aobasis.o: ../../src/libxtp/aobasis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aobasis.o ../../src/libxtp/aobasis.cc

${OBJECTDIR}/_ext/2a3bfc3d/aocoulomb_g.o: ../../src/libxtp/aomatrices/aocoulomb_g.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aocoulomb_g.o ../../src/libxtp/aomatrices/aocoulomb_g.cc

${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o: ../../src/libxtp/aomatrices/aodipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o ../../src/libxtp/aomatrices/aodipole.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o: ../../src/libxtp/aomatrices/aoecp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o ../../src/libxtp/aomatrices/aoecp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o: ../../src/libxtp/aomatrices/aoesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o ../../src/libxtp/aomatrices/aoesp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o: ../../src/libxtp/aomatrices/aokinetic.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o ../../src/libxtp/aomatrices/aokinetic.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o: ../../src/libxtp/aomatrices/aomatrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o ../../src/libxtp/aomatrices/aomatrix.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o: ../../src/libxtp/aomatrices/aomomentum.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o ../../src/libxtp/aomatrices/aomomentum.cc

${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o: ../../src/libxtp/aomatrices/aooverlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o ../../src/libxtp/aomatrices/aooverlap.cc

${OBJECTDIR}/_ext/1ce08c7a/aoshell.o: ../../src/libxtp/aoshell.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aoshell.o ../../src/libxtp/aoshell.cc

${OBJECTDIR}/_ext/1ce08c7a/apolarsite.o: ../../src/libxtp/apolarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/apolarsite.o ../../src/libxtp/apolarsite.cc

${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o: ../../src/libxtp/bsecoupling.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o ../../src/libxtp/bsecoupling.cc

${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o: ../../src/libxtp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o ../../src/libxtp/calculatorfactory.cc

${OBJECTDIR}/_ext/50c6f79c/eoutersphere.o: ../../src/libxtp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/50c6f79c
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/50c6f79c/eoutersphere.o ../../src/libxtp/calculators/eoutersphere.cc

${OBJECTDIR}/_ext/50c6f79c/jobwriter.o: ../../src/libxtp/calculators/jobwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/50c6f79c
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/50c6f79c/jobwriter.o ../../src/libxtp/calculators/jobwriter.cc

${OBJECTDIR}/_ext/5e78919f/dftengine.o: ../../src/libxtp/dftengine/dftengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/5e78919f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/5e78919f/dftengine.o ../../src/libxtp/dftengine/dftengine.cc

${OBJECTDIR}/_ext/1ce08c7a/espfit.o: ../../src/libxtp/espfit.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/espfit.o ../../src/libxtp/espfit.cc

${OBJECTDIR}/_ext/1ce08c7a/ewaldactor.o: ../../src/libxtp/ewaldactor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/ewaldactor.o ../../src/libxtp/ewaldactor.cc

${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o: ../../src/libxtp/extractorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o ../../src/libxtp/extractorfactory.cc

${OBJECTDIR}/_ext/1ce08c7a/fragment.o: ../../src/libxtp/fragment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/fragment.o ../../src/libxtp/fragment.cc

${OBJECTDIR}/_ext/1ce08c7a/gdma.o: ../../src/libxtp/gdma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gdma.o ../../src/libxtp/gdma.cc

${OBJECTDIR}/_ext/1ce08c7a/grid.o: ../../src/libxtp/grid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/grid.o ../../src/libxtp/grid.cc

${OBJECTDIR}/_ext/1ce08c7a/gsl_boost_ublas_matrix_prod.o: ../../src/libxtp/gsl_boost_ublas_matrix_prod.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gsl_boost_ublas_matrix_prod.o ../../src/libxtp/gsl_boost_ublas_matrix_prod.cc

${OBJECTDIR}/_ext/6916996f/bse.o: ../../src/libxtp/gwbse/bse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/bse.o ../../src/libxtp/gwbse/bse.cc

${OBJECTDIR}/_ext/6916996f/gwa.o: ../../src/libxtp/gwbse/gwa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwa.o ../../src/libxtp/gwbse/gwa.cc

${OBJECTDIR}/_ext/6916996f/gwbse.o: ../../src/libxtp/gwbse/gwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwbse.o ../../src/libxtp/gwbse/gwbse.cc

${OBJECTDIR}/_ext/6916996f/rpa.o: ../../src/libxtp/gwbse/rpa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/rpa.o ../../src/libxtp/gwbse/rpa.cc

${OBJECTDIR}/_ext/1ce08c7a/job.o: ../../src/libxtp/job.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/job.o ../../src/libxtp/job.cc

${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o: ../../src/libxtp/jobapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o ../../src/libxtp/jobapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o: ../../src/libxtp/jobcalculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o ../../src/libxtp/jobcalculatorfactory.cc

${OBJECTDIR}/_ext/f632c409/egwbse.o: ../../src/libxtp/jobcalculators/egwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/egwbse.o ../../src/libxtp/jobcalculators/egwbse.cc

${OBJECTDIR}/_ext/f632c409/idft.o: ../../src/libxtp/jobcalculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/idft.o ../../src/libxtp/jobcalculators/idft.cc

${OBJECTDIR}/_ext/f632c409/iexcitoncl.o: ../../src/libxtp/jobcalculators/iexcitoncl.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/iexcitoncl.o ../../src/libxtp/jobcalculators/iexcitoncl.cc

${OBJECTDIR}/_ext/f632c409/igwbse.o: ../../src/libxtp/jobcalculators/igwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/igwbse.o ../../src/libxtp/jobcalculators/igwbse.cc

${OBJECTDIR}/_ext/1ce08c7a/molecule.o: ../../src/libxtp/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/molecule.o ../../src/libxtp/molecule.cc

${OBJECTDIR}/_ext/4d261038/numerical_integrations.o: ../../src/libxtp/numerical_integration/numerical_integrations.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/numerical_integrations.o ../../src/libxtp/numerical_integration/numerical_integrations.cc

${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o: ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc

${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o: ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc

${OBJECTDIR}/_ext/1ce08c7a/orbitals.o: ../../src/libxtp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/orbitals.o ../../src/libxtp/orbitals.cc

${OBJECTDIR}/_ext/1ce08c7a/overlap.o: ../../src/libxtp/overlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/overlap.o ../../src/libxtp/overlap.cc

${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o: ../../src/libxtp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/parallelpaircalc.o ../../src/libxtp/parallelpaircalc.cc

${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o: ../../src/libxtp/parallelxjobcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/parallelxjobcalc.o ../../src/libxtp/parallelxjobcalc.cc

${OBJECTDIR}/_ext/1ce08c7a/poissongrid.o: ../../src/libxtp/poissongrid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/poissongrid.o ../../src/libxtp/poissongrid.cc

${OBJECTDIR}/_ext/1ce08c7a/polarbackground.o: ../../src/libxtp/polarbackground.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/polarbackground.o ../../src/libxtp/polarbackground.cc

${OBJECTDIR}/_ext/1ce08c7a/polarfrag.o: ../../src/libxtp/polarfrag.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/polarfrag.o ../../src/libxtp/polarfrag.cc

${OBJECTDIR}/_ext/1ce08c7a/polarseg.o: ../../src/libxtp/polarseg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/polarseg.o ../../src/libxtp/polarseg.cc

${OBJECTDIR}/_ext/1ce08c7a/polarsite.o: ../../src/libxtp/polarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/polarsite.o ../../src/libxtp/polarsite.cc

${OBJECTDIR}/_ext/1ce08c7a/polartop.o: ../../src/libxtp/polartop.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/polartop.o ../../src/libxtp/polartop.cc

${OBJECTDIR}/_ext/1ce08c7a/progressobserver.o: ../../src/libxtp/progressobserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/progressobserver.o ../../src/libxtp/progressobserver.cc

${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o: ../../src/libxtp/qmapemachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o ../../src/libxtp/qmapemachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o: ../../src/libxtp/qmcalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmcalculator.o ../../src/libxtp/qmcalculator.cc

${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o: ../../src/libxtp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o ../../src/libxtp/qmdatabase.cc

${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o: ../../src/libxtp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o ../../src/libxtp/qmmachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmnblist.o: ../../src/libxtp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmnblist.o ../../src/libxtp/qmnblist.cc

${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o: ../../src/libxtp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o ../../src/libxtp/qmpackagefactory.cc

${OBJECTDIR}/_ext/60851bbe/gaussian.o: ../../src/libxtp/qmpackages/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/gaussian.o ../../src/libxtp/qmpackages/gaussian.cc

${OBJECTDIR}/_ext/60851bbe/nwchem.o: ../../src/libxtp/qmpackages/nwchem.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/nwchem.o ../../src/libxtp/qmpackages/nwchem.cc

${OBJECTDIR}/_ext/60851bbe/orca.o: ../../src/libxtp/qmpackages/orca.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/orca.o ../../src/libxtp/qmpackages/orca.cc

${OBJECTDIR}/_ext/60851bbe/turbomole.o: ../../src/libxtp/qmpackages/turbomole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/turbomole.o ../../src/libxtp/qmpackages/turbomole.cc

${OBJECTDIR}/_ext/1ce08c7a/qmpair.o: ../../src/libxtp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmpair.o ../../src/libxtp/qmpair.cc

${OBJECTDIR}/_ext/1ce08c7a/qmtool.o: ../../src/libxtp/qmtool.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmtool.o ../../src/libxtp/qmtool.cc

${OBJECTDIR}/_ext/1ce08c7a/segment.o: ../../src/libxtp/segment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/segment.o ../../src/libxtp/segment.cc

${OBJECTDIR}/_ext/1ce08c7a/segmenttype.o: ../../src/libxtp/segmenttype.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/segmenttype.o ../../src/libxtp/segmenttype.cc

${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o: ../../src/libxtp/sqlapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o ../../src/libxtp/sqlapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o: ../../src/libxtp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o ../../src/libxtp/statesaversqlite.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o: ../../src/libxtp/threecenter_rep.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o ../../src/libxtp/threecenter_rep.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters.o: ../../src/libxtp/threecenters.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters.o ../../src/libxtp/threecenters.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o: ../../src/libxtp/threecenters_dft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o ../../src/libxtp/threecenters_dft.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o: ../../src/libxtp/threecenters_g.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_g.o ../../src/libxtp/threecenters_g.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o: ../../src/libxtp/threecenters_tools.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o ../../src/libxtp/threecenters_tools.cc

${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o: ../../src/libxtp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o ../../src/libxtp/toolfactory.cc

${OBJECTDIR}/_ext/69ca5806/molpol.o: ../../src/libxtp/tools/molpol.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/69ca5806
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/69ca5806/molpol.o ../../src/libxtp/tools/molpol.cc

${OBJECTDIR}/_ext/1ce08c7a/topology.o: ../../src/libxtp/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/topology.o ../../src/libxtp/topology.cc

${OBJECTDIR}/_ext/1ce08c7a/version.o: ../../src/libxtp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version.o ../../src/libxtp/version.cc

${OBJECTDIR}/_ext/1ce08c7a/version_nb.o: ../../src/libxtp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version_nb.o ../../src/libxtp/version_nb.cc

${OBJECTDIR}/_ext/1ce08c7a/xinductor.o: ../../src/libxtp/xinductor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xinductor.o ../../src/libxtp/xinductor.cc

${OBJECTDIR}/_ext/1ce08c7a/xinteractor.o: ../../src/libxtp/xinteractor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xinteractor.o ../../src/libxtp/xinteractor.cc

${OBJECTDIR}/_ext/1ce08c7a/xjob.o: ../../src/libxtp/xjob.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xjob.o ../../src/libxtp/xjob.cc

${OBJECTDIR}/_ext/1ce08c7a/xmapper.o: ../../src/libxtp/xmapper.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xmapper.o ../../src/libxtp/xmapper.cc

${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o: ../../src/libxtp/xtpapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o ../../src/libxtp/xtpapplication.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ../../src/libxtp/libxtp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
