/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <boost/format.hpp>
#include <votca/tools/constants.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/bsecoupling.h>

namespace votca {
namespace xtp {

using namespace std;
using boost::format;
using namespace tools;

void BSECoupling::Initialize(Property& options) {

#if (GWBSE_DOUBLE)
  XTP_LOG(logDEBUG, *_pLog) << " Compiled with full double support" << flush;
#else
  XTP_LOG(logDEBUG, *_pLog)
      << " Compiled with float/double mixture (standard)" << flush;
#endif

  std::string key = Identify();
  _doSinglets = false;
  _doTriplets = false;
  _output_perturbation = false;

  _openmp_threads = 0;

  if (options.exists(key + ".openmp")) {
    _openmp_threads = options.get(key + ".openmp").as<int>();
  }

  string spintype = options.get(key + ".spin").as<string>();
  if (spintype == "all") {
    _doSinglets = true;
    _doTriplets = true;
  } else if (spintype == "triplet") {
    _doTriplets = true;
  } else if (spintype == "singlet") {
    _doSinglets = true;
  } else {
    throw std::runtime_error(
        (boost::format(
             "Choice % for type not known. Available singlet,triplet,all") %
         spintype)
            .str());
  }

  if (options.exists(key + ".algorithm")) {
    string algorithm = options.get(key + ".algorithm").as<string>();
    if (algorithm == "perturbation") {
      _output_perturbation = true;
    }
  }

  _levA = options.get(key + ".moleculeA.states").as<int>();
  _levB = options.get(key + ".moleculeB.states").as<int>();
  _occA = options.get(key + ".moleculeA.occLevels").as<int>();
  _occB = options.get(key + ".moleculeB.occLevels").as<int>();
  _unoccA = options.get(key + ".moleculeA.unoccLevels").as<int>();
  _unoccB = options.get(key + ".moleculeB.unoccLevels").as<int>();
}

void BSECoupling::WriteToProperty(const Orbitals& orbitalsA,
                                  const Orbitals& orbitalsB, Property& summary,
                                  const QMState& stateA,
                                  const QMState& stateB) {
  Property& coupling_summary = summary.add("coupling", "");
  double energyA = 0;
  double energyB = 0;
  double JAB_pert = 0;
  double JAB_diag = 0;
  if (stateA.Type() == QMStateType::Singlet) {
    energyA = orbitalsA.BSESingletEnergies()(stateA.Index()) * conv::hrt2ev;
    energyB = orbitalsB.BSESingletEnergies()(stateB.Index()) * conv::hrt2ev;
    JAB_pert = getSingletCouplingElement(stateA.Index(), stateB.Index(), 0);
    JAB_diag = getSingletCouplingElement(stateA.Index(), stateB.Index(), 1);
  } else if (stateA.Type() == QMStateType::Triplet) {
    energyA = orbitalsA.BSETripletEnergies()(stateA.Index()) * conv::hrt2ev;
    energyB = orbitalsB.BSETripletEnergies()(stateB.Index()) * conv::hrt2ev;
    JAB_pert = getTripletCouplingElement(stateA.Index(), stateB.Index(), 0);
    JAB_diag = getTripletCouplingElement(stateA.Index(), stateB.Index(), 1);
  }
  coupling_summary.setAttribute("stateA", stateA.ToString());
  coupling_summary.setAttribute("stateB", stateB.ToString());
  coupling_summary.setAttribute("eA", (format("%1$1.6e") % energyA).str());
  coupling_summary.setAttribute("eB", (format("%1$1.6e") % energyB).str());
  coupling_summary.setAttribute("j_pert", (format("%1$1.6e") % JAB_pert).str());
  coupling_summary.setAttribute("j_diag", (format("%1$1.6e") % JAB_diag).str());
}

void BSECoupling::Addoutput(Property& type_summary, const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) {
  tools::Property& bsecoupling = type_summary.add(Identify(), "");
  string algorithm = "j_diag";
  if (_output_perturbation) {
    algorithm = "j_pert";
  }
  if (_doSinglets) {
    QMStateType singlet = QMStateType(QMStateType::Singlet);
    Property& singlet_summary = bsecoupling.add(singlet.ToLongString(), "");
    singlet_summary.setAttribute("algorithm", algorithm);
    for (int stateA = 0; stateA < _levA; ++stateA) {
      QMState qmstateA = QMState(singlet, stateA, false);
      for (int stateB = 0; stateB < _levB; ++stateB) {
        QMState qmstateB = QMState(singlet, stateB, false);
        WriteToProperty(orbitalsA, orbitalsB, singlet_summary, qmstateA,
                        qmstateB);
      }
    }
  }

  if (_doTriplets) {
    QMStateType triplet = QMStateType(QMStateType::Triplet);
    Property& triplet_summary = bsecoupling.add(triplet.ToLongString(), "");
    triplet_summary.setAttribute("algorithm", algorithm);
    for (int stateA = 0; stateA < _levA; ++stateA) {
      QMState qmstateA = QMState(triplet, stateA, false);
      for (int stateB = 0; stateB < _levB; ++stateB) {
        QMState qmstateB = QMState(triplet, stateB, false);
        WriteToProperty(orbitalsA, orbitalsB, triplet_summary, qmstateA,
                        qmstateB);
      }
    }
  }
}

double BSECoupling::getSingletCouplingElement(int levelA, int levelB,
                                              int methodindex) {
  return JAB_singlet[methodindex](levelA, levelB + _levA) *
         votca::tools::conv::hrt2ev;
}

double BSECoupling::getTripletCouplingElement(int levelA, int levelB,
                                              int methodindex) {
  return JAB_triplet[methodindex](levelA, levelB + _levA) *
         votca::tools::conv::hrt2ev;
}

/**
 * \brief evaluates electronic couplings
 *
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 */
void BSECoupling::CalculateCouplings(const Orbitals& orbitalsA,
                                     const Orbitals& orbitalsB,
                                     Orbitals& orbitalsAB) {
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "  Calculating exciton couplings" << flush;
#ifdef _OPENMP

  if (_openmp_threads > 0) omp_set_num_threads(_openmp_threads);
  XTP_LOG(logDEBUG, *_pLog) << TimeStamp() << " Using " << omp_get_max_threads()
                            << " threads" << flush;
#endif

  CheckAtomCoordinates(orbitalsA, orbitalsB, orbitalsAB);

  // constructing the direct product orbA x orbB
  int basisA = orbitalsA.getBasisSetSize();
  int basisB = orbitalsB.getBasisSetSize();

  if ((basisA == 0) || (basisB == 0)) {
    throw std::runtime_error("Basis set size is not stored in monomers");
  }

  // number of levels stored in monomers
  int levelsA = orbitalsA.getBasisSetSize();
  int levelsB = orbitalsB.getBasisSetSize();

  // get exciton information of molecule A
  int _bseA_cmax = orbitalsA.getBSEcmax();
  int _bseA_cmin = orbitalsA.getBSEcmin();
  int _bseA_vmax = orbitalsA.getBSEvmax();
  int _bseA_vmin = orbitalsA.getBSEvmin();
  int _bseA_vtotal = _bseA_vmax - _bseA_vmin + 1;
  int _bseA_ctotal = _bseA_cmax - _bseA_cmin + 1;
  int _bseA_size = _bseA_vtotal * _bseA_ctotal;
  int _bseA_singlet_exc = orbitalsA.BSESingletCoefficients().cols();
  int _bseA_triplet_exc = orbitalsA.BSETripletCoefficients().cols();

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   molecule A has " << _bseA_singlet_exc
      << " singlet excitons with dimension " << _bseA_size << flush;
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   molecule A has " << _bseA_triplet_exc
      << " triplet excitons with dimension " << _bseA_size << flush;

  // now, two storage assignment matrices for two-particle functions
  Eigen::MatrixXi combA;
  combA.resize(_bseA_size, 2);
  int cnt = 0;
  for (int _v = 0; _v < _bseA_vtotal; _v++) {
    for (int _c = 0; _c < _bseA_ctotal; _c++) {
      combA(cnt, 0) = _v;
      combA(cnt, 1) = _bseA_vtotal + _c;
      cnt++;
    }
  }

  // get exciton information of molecule B
  int _bseB_cmax = orbitalsB.getBSEcmax();
  int _bseB_cmin = orbitalsB.getBSEcmin();
  int _bseB_vmax = orbitalsB.getBSEvmax();
  int _bseB_vmin = orbitalsB.getBSEvmin();
  int _bseB_vtotal = _bseB_vmax - _bseB_vmin + 1;
  int _bseB_ctotal = _bseB_cmax - _bseB_cmin + 1;
  int _bseB_size = _bseB_vtotal * _bseB_ctotal;
  int _bseB_singlet_exc = orbitalsB.BSESingletCoefficients().cols();
  int _bseB_triplet_exc = orbitalsB.BSETripletCoefficients().cols();

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   molecule B has " << _bseB_singlet_exc
      << " singlet excitons with dimension " << _bseB_size << flush;
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   molecule B has " << _bseB_triplet_exc
      << " triplet excitons with dimension " << _bseB_size << flush;

  // now, two storage assignment matrices for two-particle functions
  Eigen::MatrixXi combB;
  combB.resize(_bseB_size, 2);
  cnt = 0;
  for (int _v = 0; _v < _bseB_vtotal; _v++) {
    for (int _c = 0; _c < _bseB_ctotal; _c++) {
      combB(cnt, 0) = _bseA_vtotal + _bseA_ctotal + _v;
      combB(cnt, 1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;
      cnt++;
    }
  }

  if (_levA > _bseA_singlet_exc) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of excitons you want is greater than stored for molecule "
           "A. Setting to max number available"
        << flush;
    _levA = _bseA_singlet_exc;
  }
  if (_levB > _bseB_singlet_exc) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of excitons you want is greater than stored for molecule "
           "B. Setting to max number available"
        << flush;
    _levB = _bseB_singlet_exc;
  }

  if (_levA > _bseA_singlet_exc) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of Frenkel states you want is greater than stored for "
           "molecule A. Setting to max number available"
        << flush;
    _levA = _bseA_singlet_exc;
  }
  if (_levB > _bseB_singlet_exc) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of Frenkel states you want is greater than stored for "
           "molecule B. Setting to max number available"
        << flush;
    _levB = _bseB_singlet_exc;
  }

  if (_unoccA > _bseA_ctotal) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of occupied orbitals in molecule A for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    _unoccA = _bseA_ctotal;
  } else if (_unoccA < 0) {
    _unoccA = _bseA_ctotal;
  }
  if (_unoccB > _bseB_ctotal) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of occupied orbitals in molecule B for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    _unoccB = _bseB_ctotal;
  } else if (_unoccB < 0) {
    _unoccB = _bseB_ctotal;
  }

  if (_occA > _bseA_vtotal) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of unoccupied orbitals in molecule A for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    _occA = _bseA_vtotal;
  } else if (_occA < 0) {
    _occA = _bseA_vtotal;
  }
  if (_occB > _bseB_vtotal) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "  Number of unoccupied orbitals in molecule B for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    _occB = _bseB_vtotal;
  } else if (_occB < 0) {
    _occB = _bseB_vtotal;
  }

  // get exciton information of pair AB
  int _bseAB_cmax = orbitalsAB.getBSEcmax();
  int _bseAB_cmin = orbitalsAB.getBSEcmin();
  int _bseAB_vmax = orbitalsAB.getBSEvmax();
  int _bseAB_vmin = orbitalsAB.getBSEvmin();
  int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin + 1;
  int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin + 1;
  int _bseAB_size = _bseAB_vtotal * _bseAB_ctotal;
  // check if electron-hole interaction matrices are stored
  if (!orbitalsAB.hasEHinteraction_triplet() && _doTriplets) {
    throw std::runtime_error("BSE EH for triplets not stored ");
  }
  if (!orbitalsAB.hasEHinteraction_singlet() && _doSinglets) {
    throw std::runtime_error("BSE EH for singlets not stored ");
  }
  const MatrixXfd& eh_t = orbitalsAB.eh_t();
  const MatrixXfd& eh_s = orbitalsAB.eh_s();
  if (_doTriplets) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "   dimer AB has BSE EH interaction triplet with dimension "
        << eh_t.rows() << " x " << eh_t.cols() << flush;
  }
  if (_doSinglets) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp()
        << "   dimer AB has BSE EH interaction singlet with dimension "
        << eh_s.rows() << " x " << eh_s.cols() << flush;
  }
  // now, two storage assignment matrices for two-particle functions
  Eigen::MatrixXi combAB;
  combAB.resize(_bseAB_size, 2);
  cnt = 0;
  for (int _v = 0; _v < _bseAB_vtotal; _v++) {
    for (int _c = 0; _c < _bseAB_ctotal; _c++) {
      combAB(cnt, 0) = _bseAB_vmin + _v;
      combAB(cnt, 1) = _bseAB_vmin + _bseAB_vtotal + _c;
      cnt++;
    }
  }

  // DFT levels of monomers can be reduced to those used in BSE
  levelsA = _bseA_vtotal + _bseA_ctotal;
  levelsB = _bseB_vtotal + _bseB_ctotal;
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   levels used in BSE of molA: " << _bseA_vmin
      << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal
      << flush;
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   levels used in BSE of molB: " << _bseB_vmin
      << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal
      << flush;

  if ((levelsA == 0) || (levelsB == 0)) {
    throw std::runtime_error(
        "No information about number of occupied/unoccupied levels is stored");
  }

  //       | Orbitals_A          0 |      | Overlap_A |
  //       | 0          Orbitals_B |.T  X   | Overlap_B |  X  ( Orbitals_AB )

  Eigen::MatrixXd psi_AxB =
      Eigen::MatrixXd::Zero(levelsA + levelsB, basisA + basisB);
  // constructing merged orbitals
  psi_AxB.block(0, 0, levelsA, basisA) = orbitalsA.MOCoefficients().block(
      _bseA_vmin, 0, _bseA_cmax + 1 - _bseA_vmin, basisA);
  psi_AxB.block(levelsA, basisA, levelsB, basisB) =
      orbitalsB.MOCoefficients().block(_bseB_vmin, 0,
                                       _bseB_cmax + 1 - _bseA_vmin, basisB);

  // psi_AxB * S_AB * psi_AB
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   projecting monomer onto dimer orbitals" << flush;

  Eigen::MatrixXd overlapAB;
  if (orbitalsAB.hasAOOverlap()) {
    XTP_LOG(logDEBUG, *_pLog)
        << "Reading overlap matrix from orbitals" << flush;
    overlapAB = orbitalsAB.AOOverlap();
  } else {
    XTP_LOG(logDEBUG, *_pLog) << "Calculating overlap matrix for basisset: "
                              << orbitalsAB.getDFTbasisName() << flush;
    overlapAB = CalculateOverlapMatrix(orbitalsAB);
  }

  Eigen::MatrixXd psi_AxB_dimer_basis =
      psi_AxB.transpose() * overlapAB * orbitalsAB.MOCoefficients();
  overlapAB.resize(0, 0);
  int LevelsA = levelsA;
  for (int i = 0; i < psi_AxB_dimer_basis.rows(); i++) {
    double mag = psi_AxB_dimer_basis.row(i).squaredNorm();
    if (mag < 0.95) {
      int monomer = 0;
      int level = 0;
      if (i < LevelsA) {
        monomer = 1;
        level = _bseA_vmin + i;
      } else {
        monomer = 2;
        level = _bseB_vmin + i - levelsA;
      }
      XTP_LOG(logERROR, *_pLog) << "\nERROR: " << i << " Projection of orbital "
                                << level << " of monomer " << monomer
                                << " on dimer is insufficient,mag=" << mag
                                << " maybe the orbital order is screwed up, "
                                   "otherwise increase dimer basis.\n"
                                << flush;
    }
  }

  // notation AB is CT states with A+B-, BA is the counterpart
  // Setting up CT-states:
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   Setting up CT-states" << flush;
  // Number of A+B- states
  int noAB = _occA * _unoccB;
  // Number of A-B+ states
  int noBA = _unoccA * _occB;

  Eigen::MatrixXi comb_CTAB = Eigen::MatrixXi::Zero(noAB, 2);
  cnt = 0;
  // iterate A over occupied, B over unoccupied
  int v_start = _bseA_vtotal - _occA;
  for (int _v = v_start; _v < _bseA_vtotal; _v++) {
    for (int _c = 0; _c < _unoccB; _c++) {
      comb_CTAB(cnt, 0) = _v;
      comb_CTAB(cnt, 1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;

      cnt++;
    }
  }
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "  " << noAB << " CT states A+B- created" << flush;

  Eigen::MatrixXi comb_CTBA = Eigen::MatrixXi::Zero(noBA, 2);
  cnt = 0;
  // iterate A over unoccupied, B over occupied
  v_start = _bseB_vtotal - _occB;
  for (int _v = v_start; _v < _bseB_vtotal; _v++) {
    for (int _c = 0; _c < _unoccA; _c++) {
      comb_CTBA(cnt, 0) = _bseA_vtotal + _bseA_ctotal + _v;
      comb_CTBA(cnt, 1) = _bseA_vtotal + _c;

      cnt++;
    }
  }
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "  " << noBA << " CT states B+A- created" << flush;

  // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the
  // i-th excitation

  ctAB.resize(noAB, _bseAB_size);
#pragma omp parallel for
  for (int _i_CT = 0; _i_CT < noAB; _i_CT++) {
    for (int _i_bseAB = 0; _i_bseAB < _bseAB_size; _i_bseAB++) {
      ctAB(_i_CT, _i_bseAB) =
          psi_AxB_dimer_basis(comb_CTAB(_i_CT, 0), combAB(_i_bseAB, 0)) *
          psi_AxB_dimer_basis(comb_CTAB(_i_CT, 1), combAB(_i_bseAB, 1));
    }
  }

  ctBA.resize(noBA, _bseAB_size);
#pragma omp parallel for
  for (int _i_CT = 0; _i_CT < noBA; _i_CT++) {
    for (int _i_bseAB = 0; _i_bseAB < _bseAB_size; _i_bseAB++) {
      ctBA(_i_CT, _i_bseAB) =
          psi_AxB_dimer_basis(comb_CTBA(_i_CT, 0), combAB(_i_bseAB, 0)) *
          psi_AxB_dimer_basis(comb_CTBA(_i_CT, 1), combAB(_i_bseAB, 1));
    }
  }

  _kap.resize(_bseA_size, _bseAB_size);
#pragma omp parallel for
  for (int _i_bseA = 0; _i_bseA < _bseA_size; _i_bseA++) {
    for (int _i_bseAB = 0; _i_bseAB < _bseAB_size; _i_bseAB++) {
      _kap(_i_bseA, _i_bseAB) =
          psi_AxB_dimer_basis(combA(_i_bseA, 0), combAB(_i_bseAB, 0)) *
          psi_AxB_dimer_basis(combA(_i_bseA, 1), combAB(_i_bseAB, 1));
    }
  }

  _kbp.resize(_bseB_size, _bseAB_size);
#pragma omp parallel for
  for (int _i_bseB = 0; _i_bseB < _bseB_size; _i_bseB++) {
    for (int _i_bseAB = 0; _i_bseAB < _bseAB_size; _i_bseAB++) {
      _kbp(_i_bseB, _i_bseAB) =
          psi_AxB_dimer_basis(combB(_i_bseB, 0), combAB(_i_bseAB, 0)) *
          psi_AxB_dimer_basis(combB(_i_bseB, 1), combAB(_i_bseAB, 1));
    }
  }

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   construct projection of product functions "
      << flush;

  psi_AxB_dimer_basis.resize(0, 0);
  combAB.resize(0, 0);
  combA.resize(0, 0);
  combB.resize(0, 0);
  // now the different spin types
  if (_doSinglets) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << "   Evaluating singlets" << flush;
    Eigen::MatrixXd Hamiltonian_AB = eh_s.cast<double>();
    XTP_LOG(logDEBUG, *_pLog) << TimeStamp() << "   Setup Hamiltonian" << flush;
    const Eigen::MatrixXd bseA_T =
        orbitalsA.BSESingletCoefficients()
            .block(0, 0, orbitalsA.BSESingletCoefficients().rows(), _levA)
            .transpose()
            .cast<double>();
    const Eigen::MatrixXd bseB_T =
        orbitalsB.BSESingletCoefficients()
            .block(0, 0, orbitalsB.BSESingletCoefficients().rows(), _levB)
            .transpose()
            .cast<double>();

    JAB_singlet = ProjectExcitons(bseA_T, bseB_T, Hamiltonian_AB);
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << "   calculated singlet couplings " << flush;
  }

  if (_doTriplets) {
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << "   Evaluating triplets" << flush;
    Eigen::MatrixXd Hamiltonian_AB = eh_s.cast<double>();
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << "  Converted Hamiltonian to double" << flush;
    const Eigen::MatrixXd bseA_T =
        orbitalsA.BSETripletCoefficients()
            .block(0, 0, orbitalsA.BSETripletCoefficients().rows(), _levA)
            .transpose()
            .cast<double>();
    const Eigen::MatrixXd bseB_T =
        orbitalsB.BSETripletCoefficients()
            .block(0, 0, orbitalsB.BSETripletCoefficients().rows(), _levB)
            .transpose()
            .cast<double>();
    JAB_triplet = ProjectExcitons(bseA_T, bseB_T, Hamiltonian_AB);
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << "   calculated triplet couplings " << flush;
  }

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "  Done with exciton couplings" << flush;
  return;
};

std::vector<Eigen::MatrixXd> BSECoupling::ProjectExcitons(
    const Eigen::MatrixXd& bseA_T, const Eigen::MatrixXd& bseB_T,
    Eigen::MatrixXd& H) {

  // get projection of monomer excitons on dimer product functions
  Eigen::MatrixXd _proj_excA = bseA_T * _kap;
  Eigen::MatrixXd _proj_excB = bseB_T * _kbp;

  _bse_exc = _levA + _levB;
  int ctABsize = ctAB.rows();
  int ctBAsize = ctBA.rows();
  _ct = ctABsize + ctBAsize;
  int nobasisfunc = H.rows();

  Eigen::MatrixXd fe_states = Eigen::MatrixXd::Zero(_bse_exc, nobasisfunc);
  fe_states.block(0, 0, _levA, nobasisfunc) = _proj_excA;
  fe_states.block(_levA, 0, _levB, nobasisfunc) = _proj_excB;

  Eigen::MatrixXd ct_states = Eigen::MatrixXd::Zero(_ct, nobasisfunc);
  if (_ct > 0) {
    // orthogonalize ct-states with respect to the FE states.
    XTP_LOG(logDEBUG, *_pLog)
        << TimeStamp() << " Orthogonalizing CT-states with respect to FE-states"
        << flush;

    if (ctABsize > 0) {
      ct_states.block(0, 0, ctABsize, nobasisfunc) = ctAB;
    }
    if (ctBAsize > 0) {
      ct_states.block(ctABsize, 0, ctBAsize, nobasisfunc) = ctBA;
    }

    // orthogonalize ct-states with respect to FE states
    Eigen::MatrixXd correction = ct_states * fe_states.transpose() * fe_states;
    ct_states -= correction;
    correction.resize(0, 0);
    // normalize
    Eigen::VectorXd norm = ct_states.rowwise().norm();
    for (int i = 0; i < _ct; i++) {
      ct_states.row(i) /= norm(i);
    }
    int minstateindex = 0;
    double minnorm = norm.minCoeff(&minstateindex);
    if (minnorm < 0.95) {
      XTP_LOG(logDEBUG, *_pLog)
          << TimeStamp() << " WARNING: CT-state " << minstateindex
          << " norm is only " << minnorm << flush;
    }
  }
  Eigen::MatrixXd projection(_bse_exc + _ct, nobasisfunc);
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << " merging projections into one vector  " << flush;
  projection.block(0, 0, _bse_exc, nobasisfunc) = fe_states;

  if (_ct > 0) {
    projection.block(_bse_exc, 0, _ct, nobasisfunc) = ct_states;
  }
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   Setting up coupling matrix size " << _bse_exc + _ct
      << "x" << _bse_exc + _ct << flush;
  // matrix _J
  //  E_A         J_AB        J_A_ABCT        J_A_BACT
  //  J_BA        E_B         J_B_ABCT        J_B_BACT
  //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
  //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT

  // this only works for hermitian/symmetric H so only in TDA

  Eigen::MatrixXd J_dimer = projection * H * projection.transpose();
  H.resize(0, 0);

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   Setting up overlap matrix size " << _bse_exc + _ct
      << "x" << _bse_exc + _ct << flush;
  Eigen::MatrixXd S_dimer = projection * projection.transpose();

  projection.resize(0, 0);
  if (tools::globals::verbose && _bse_exc + _ct < 100) {
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
    XTP_LOG(logDEBUG, *_pLog) << "_J_dimer[Ryd]" << flush;

    XTP_LOG(logDEBUG, *_pLog) << J_dimer << flush;
    XTP_LOG(logDEBUG, *_pLog) << "_S_dimer" << flush;

    XTP_LOG(logDEBUG, *_pLog) << S_dimer << flush;
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_dimer);
  Eigen::MatrixXd Sm1 = es.operatorInverseSqrt();
  J_dimer = Sm1 * J_dimer * Sm1;

  if (tools::globals::verbose && _bse_exc + _ct < 100) {
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
    XTP_LOG(logDEBUG, *_pLog) << "_J_ortho[Ryd]" << flush;
    XTP_LOG(logDEBUG, *_pLog) << J_dimer << flush;
    XTP_LOG(logDEBUG, *_pLog) << "_S-1/2" << flush;
    XTP_LOG(logDEBUG, *_pLog) << Sm1 << flush;
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
  }
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   Smallest value of dimer overlapmatrix is "
      << es.eigenvalues()(0) << flush;

  std::vector<Eigen::MatrixXd> J;

  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "   Running Perturbation algorithm" << flush;
  J.push_back(Perturbation(J_dimer));
  XTP_LOG(logDEBUG, *_pLog)
      << TimeStamp() << "    Running Projection algorithm" << flush;
  J.push_back(Fulldiag(J_dimer));

  if (tools::globals::verbose) {
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
    XTP_LOG(logDEBUG, *_pLog) << "Jeff_pert[Hrt]" << flush;
    XTP_LOG(logDEBUG, *_pLog) << J[0] << flush;
    XTP_LOG(logDEBUG, *_pLog) << "Jeff_diag[Hrt]" << flush;
    XTP_LOG(logDEBUG, *_pLog) << J[1] << flush;
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
  }

  return J;
}

Eigen::MatrixXd BSECoupling::Perturbation(const Eigen::MatrixXd& J_dimer) {

  Eigen::MatrixXd Jmatrix = Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);
  bool diag_ct = true;
  Eigen::MatrixXd J_result = J_dimer;
  if (_ct > 0 && diag_ct) {
    Eigen::MatrixXd transformation =
        Eigen::MatrixXd::Identity(_bse_exc + _ct, _bse_exc + _ct);
    Eigen::MatrixXd Ct = J_dimer.block(_bse_exc, _bse_exc, _ct, _ct);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ct);
    transformation.block(_bse_exc, _bse_exc, _ct, _ct) = es.eigenvectors();
    Ct.resize(0, 0);

    if (tools::globals::verbose) {
      XTP_LOG(logDEBUG, *_pLog) << "FE state hamiltonian" << flush;
      XTP_LOG(logDEBUG, *_pLog)
          << J_dimer.block(0, 0, _bse_exc, _bse_exc) << flush;
      if (_ct > 0) {
        XTP_LOG(logDEBUG, *_pLog) << "eigenvalues of CT states" << flush;
        XTP_LOG(logDEBUG, *_pLog) << es.eigenvalues() << flush;
      }
    }

    J_result = transformation.transpose() * J_dimer * transformation;
    if (tools::globals::verbose && _bse_exc + _ct < 100) {
      XTP_LOG(logDEBUG, *_pLog)
          << "---------------------------------------" << flush;
      XTP_LOG(logDEBUG, *_pLog) << "_J_ortho[Hrt] CT-state diag" << flush;
      XTP_LOG(logDEBUG, *_pLog) << J_result << flush;
      XTP_LOG(logDEBUG, *_pLog)
          << "---------------------------------------" << flush;
    }
  }
  for (int stateA = 0; stateA < _levA; stateA++) {
    double Ea = J_result(stateA, stateA);
    for (int stateB = 0; stateB < _levB; stateB++) {
      int stateBd = stateB + _levA;
      XTP_LOG(logDEBUG, *_pLog)
          << TimeStamp() << "   Calculating coupling between exciton A"
          << stateA + 1 << " and exciton B" << stateB + 1 << flush;
      double J = J_result(stateA, stateBd);

      double Eb = J_result(stateBd, stateBd);
      for (int k = _bse_exc; k < (_bse_exc + _ct); k++) {
        double Eab = J_result(k, k);
        if (std::abs(Eab - Ea) < 0.001) {
          XTP_LOG(logDEBUG, *_pLog)
              << TimeStamp() << "Energydifference between state A "
              << stateA + 1 << "and CT state " << k + 1 << " is " << Eab - Ea
              << "[Hrt]" << flush;
        }
        if (std::abs(Eab - Eb) < 0.001) {
          XTP_LOG(logDEBUG, *_pLog)
              << TimeStamp() << "Energydifference between state B "
              << stateB + 1 << "and CT state " << k + 1 << " is " << Eab - Eb
              << "[Hrt]" << flush;
        }
        J += 0.5 * J_result(k, stateA) * J_result(k, stateBd) *
             (1 / (Ea - Eab) + 1 / (Eb - Eab));  // Have no clue why 0.5
      }
      Jmatrix(stateA, stateBd) = J;
      Jmatrix(stateBd, stateA) = J;
    }
  }
  return Jmatrix;
}

Eigen::MatrixXd BSECoupling::Fulldiag(const Eigen::MatrixXd& J_dimer) {
  Eigen::MatrixXd Jmat = Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J_dimer);
  if (tools::globals::verbose && _bse_exc + _ct < 10) {
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
    XTP_LOG(logDEBUG, *_pLog) << "Eigenvectors of J" << flush;
    XTP_LOG(logDEBUG, *_pLog) << es.eigenvectors() << flush;
    XTP_LOG(logDEBUG, *_pLog) << "J_eigenvalues[Hrt]" << flush;
    XTP_LOG(logDEBUG, *_pLog) << es.eigenvalues() << flush;
    XTP_LOG(logDEBUG, *_pLog)
        << "---------------------------------------" << flush;
  }
  // Calculate projection on subspace for every pair of excitons separately
  for (int stateA = 0; stateA < _levA; stateA++) {
    for (int stateB = 0; stateB < _levB; stateB++) {
      int stateBd = stateB + _levA;
      XTP_LOG(logDEBUG, *_pLog)
          << TimeStamp() << "   Calculating coupling between exciton A"
          << stateA + 1 << " and exciton B" << stateB + 1 << flush;
      std::vector<int> index;
      std::vector<int> signvec;
      for (int i = 0; i < _bse_exc + _ct; i++) {
        if (i == int(stateA) || i == int(stateBd)) {
          double close = 0.0;
          int ind = 0;
          int sign = 0;
          // row
          for (int j = 0; j < _bse_exc + _ct; j++) {
            bool check = true;
            // if index i is already in index
            // should not happen but if one vector was similar to two others.
            for (unsigned l = 0; l < index.size(); l++) {
              if (j == index[l]) {
                check = false;
                break;
              }
            }
            if (check && std::abs(es.eigenvalues()(i, j)) > close) {
              ind = j;
              close = std::abs(es.eigenvalues()(i, j));
              if (es.eigenvalues()(i, j) >= 0) {
                sign = 1;
              } else {
                sign = -1;
              }
            }
          }
          index.push_back(ind);
          signvec.push_back(sign);
        }
      }

      XTP_LOG(logDEBUG, *_pLog)
          << TimeStamp() << "   Order is: [Initial state n->nth eigenvalue]"
          << flush;
      XTP_LOG(logDEBUG, *_pLog) << "    A" << stateA + 1 << ":" << stateA + 1
                                << "->" << index[0] + 1 << " ";
      XTP_LOG(logDEBUG, *_pLog) << "    B" << stateB + 1 << ":" << stateBd + 1
                                << "->" << index[1] + 1 << " " << flush;

      // setting up transformation matrix Tmat and diagonal matrix Emat for the
      // eigenvalues;
      Eigen::MatrixXd Emat = Eigen::MatrixXd::Zero(2, 2);
      Eigen::MatrixXd Tmat = Eigen::MatrixXd::Zero(2, 2);
      // find the eigenvectors which are most similar to the initial states
      // row
      for (int i = 0; i < 2; i++) {
        int k = index[i];
        double sign = signvec[i];
        double normr = 1 / std::sqrt(es.eigenvectors()(stateA, k) *
                                         es.eigenvectors()(stateA, k) +
                                     es.eigenvectors()(stateBd, k) *
                                         es.eigenvectors()(stateBd, k));
        Tmat(0, i) = sign * es.eigenvectors()(stateA, k) * normr;
        Tmat(1, i) = sign * es.eigenvectors()(stateBd, k) * normr;
        Emat(i, i) = es.eigenvectors()(k);
      }

      if ((Tmat(1, 1) * Tmat(0, 0) - Tmat(1, 0) * Tmat(0, 1)) < 0) {
        XTP_LOG(logDEBUG, *_pLog)
            << " Reduced state matrix is not in a right handed basis, "
               "multiplying second eigenvector by -1 "
            << flush;
        Tmat(0, 1) = -Tmat(0, 1);
        Tmat(1, 1) = -Tmat(1, 1);
      }

      if (tools::globals::verbose) {
        XTP_LOG(logDEBUG, *_pLog)
            << "---------------------------------------" << flush;
        XTP_LOG(logDEBUG, *_pLog) << "_T" << flush;
        XTP_LOG(logDEBUG, *_pLog) << Tmat << flush;
      }

      Eigen::MatrixXd S_small = Tmat * Tmat.transpose();
      if (tools::globals::verbose) {

        XTP_LOG(logDEBUG, *_pLog) << "S_small" << flush;
        XTP_LOG(logDEBUG, *_pLog) << S_small << flush;
      }
      // orthogonalize that matrix

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ss(S_small);
      Eigen::MatrixXd sm1 = ss.operatorInverseSqrt();
      Emat = sm1 * Emat * sm1;

      XTP_LOG(logDEBUG, *_pLog)
          << TimeStamp() << "   Smallest value of dimer overlapmatrix is "
          << ss.eigenvalues()(0) << flush;
      if (tools::globals::verbose) {

        XTP_LOG(logDEBUG, *_pLog) << "S-1/2" << flush;
        XTP_LOG(logDEBUG, *_pLog) << sm1 << flush;
        XTP_LOG(logDEBUG, *_pLog) << "E_ortho" << flush;
        XTP_LOG(logDEBUG, *_pLog) << Emat << flush;
      }
      Tmat = Tmat * sm1;

      if (tools::globals::verbose) {

        XTP_LOG(logDEBUG, *_pLog) << "T_ortho" << flush;
        XTP_LOG(logDEBUG, *_pLog) << Tmat << flush;
        XTP_LOG(logDEBUG, *_pLog)
            << "---------------------------------------" << flush;
      }

      Eigen::MatrixXd J_small = Tmat * Emat * Tmat.transpose();
      if (tools::globals::verbose) {
        XTP_LOG(logDEBUG, *_pLog) << "T_ortho*E_ortho*T_ortho^T" << flush;
        XTP_LOG(logDEBUG, *_pLog) << J_small << flush;
      }

      Jmat(stateA, stateBd) = J_small(0, 1);
      Jmat(stateBd, stateA) = J_small(1, 0);
    }
  }

  return Jmat;
}

}  // namespace xtp
}  // namespace votca
