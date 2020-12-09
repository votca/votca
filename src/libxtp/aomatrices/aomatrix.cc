/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include "votca/xtp/make_libint_work.h"
#include <vector>
// Local VOTCA includes
#include "votca/xtp/aomatrix.h"

// include libint last otherwise it overrides eigen
#include <libint2.hpp>
namespace votca {
namespace xtp {

template <libint2::Operator obtype,
          typename OperatorParams =
              typename libint2::operator_traits<obtype>::oper_params_type>
std::array<AOMatrix::MatrixLibInt, libint2::operator_traits<obtype>::nopers>
    computeOneBodyIntegrals(const AOBasis& aobasis,
                            OperatorParams oparams = OperatorParams()) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();

  std::vector<std::vector<Index>> shellpair_list = aobasis.ComputeShellPairs();

  Index nopers = static_cast<Index>(libint2::operator_traits<obtype>::nopers);
  std::array<AOMatrix::MatrixLibInt, libint2::operator_traits<obtype>::nopers>
      result;
  for (AOMatrix::MatrixLibInt& r : result) {
    r = AOMatrix::MatrixLibInt::Zero(aobasis.AOBasisSize(),
                                     aobasis.AOBasisSize());
  }

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(obtype, aobasis.getMaxNprim(),
                               static_cast<int>(aobasis.getMaxL()), 0);
  engines[0].set_params(oparams);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();

    for (Index s2 : shellpair_list[s1]) {

      engine.compute(shells[s1], shells[s2]);
      if (buf[0] == nullptr) {
        continue;  // if all integrals screened out, skip to next shell set
      }
      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();
      for (unsigned int op = 0; op != nopers; ++op) {
        Eigen::Map<const AOMatrix::MatrixLibInt> buf_mat(buf[op], n1, n2);
        result[op].block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                         // {s2,s1} block, note the transpose!
          result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
      }
    }
  }
  return result;
}

/***********************************
 * KINETIC
 ***********************************/
void AOKinetic::Fill(const AOBasis& aobasis) {
  _aomatrix = computeOneBodyIntegrals<libint2::Operator::kinetic>(aobasis)[0];
}

/***********************************
 * OVERLAP
 ***********************************/
void AOOverlap::Fill(const AOBasis& aobasis) {
  _aomatrix = computeOneBodyIntegrals<libint2::Operator::overlap>(aobasis)[0];
}

Eigen::MatrixXd AOOverlap::singleShellOverlap(const AOShell& shell) const {
  libint2::Shell::do_enforce_unit_normalization(false);
  libint2::Operator obtype = libint2::Operator::overlap;
  try {
    libint2::Engine engine(obtype, shell.getSize(),
                           static_cast<int>(shell.getL()), 0);

    const libint2::Engine::target_ptr_vec& buf = engine.results();

    libint2::Shell s = shell.LibintShell();
    engine.compute(s, s);

    Eigen::Map<const AOMatrix::MatrixLibInt> buf_mat(buf[0], shell.getNumFunc(),
                                                     shell.getNumFunc());
    return buf_mat;
  } catch (const libint2::Engine::lmax_exceeded& error) {
    std::ostringstream oss;
    oss << "\nA libint error occured:\n" << error.what() << "\n"
        << "\nYou requested a computation for a shell with angular momentum "
        << error.lmax_requested()
        << ",\nbut your libint package only supports angular momenta upto "
        << error.lmax_limit() - 1 << ".\n"
        << "\nTo fix this error you will need to reinstall libint with "
           "support\n"
           "for higher angular momenta. If you installed your own libint it\n"
           "should be reconfigured and installed with the option "
           "--with-max-am=<maxAngularMomentum>.\n"
           "If you installed libint with your OS package manager, you will\n"
           "need to setup you own libint installation with the \n"
           "--with-max-am=<maxAngularMomentum> option set."
        << std::endl;
    throw std::runtime_error(oss.str());
  }
}

Eigen::MatrixXd AOOverlap::Pseudo_InvSqrt(double etol) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());
  removedfunctions = 0;
  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  return es.eigenvectors() * diagonal.asDiagonal() *
         es.eigenvectors().transpose();
}

Eigen::MatrixXd AOOverlap::Sqrt() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  return es.operatorSqrt();
}

/***********************************
 * COULOMB
 ***********************************/
void AOCoulomb::Fill(const AOBasis& aobasis) {
  computeCoulombIntegrals(aobasis);
}

void AOCoulomb::computeCoulombIntegrals(const AOBasis& aobasis) {
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

  // build engines for each thread
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] =
      libint2::Engine(libint2::Operator::coulomb, aobasis.getMaxNprim(),
                      static_cast<int>(aobasis.getMaxL()), 0);
  engines[0].set(libint2::BraKet::xs_xs);
  for (Index i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();
    // cannot use shellpairs because this is still a two-center integral and
    // overlap screening would give wrong result
    for (Index s2 = 0; s2 <= s1; ++s2) {

      engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 0>(
          shells[s1], libint2::Shell::unit(), shells[s2],
          libint2::Shell::unit());

      if (buf[0] == nullptr) {
        continue;  // if all integrals screened out, skip to next shell set
      }
      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();

      Eigen::Map<const AOMatrix::MatrixLibInt> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                       // {s2,s1} block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  }
  _aomatrix = result;
}

// This converts V into ((S-1/2 V S-1/2)-1/2 S-1/2)T, which is needed to
// construct 4c integrals,
Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                                double etol) {

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(auxoverlap.Matrix());
  removedfunctions = 0;
  Eigen::VectorXd diagonal_overlap =
      Eigen::VectorXd::Zero(eo.eigenvalues().size());
  for (Index i = 0; i < diagonal_overlap.size(); ++i) {
    if (eo.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal_overlap(i) = 1.0 / std::sqrt(eo.eigenvalues()(i));
    }
  }
  Eigen::MatrixXd Ssqrt = eo.eigenvectors() * diagonal_overlap.asDiagonal() *
                          eo.eigenvectors().transpose();

  Eigen::MatrixXd ortho = Ssqrt * _aomatrix * Ssqrt;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());

  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  Eigen::MatrixXd Vm1 =
      es.eigenvectors() * diagonal.asDiagonal() * es.eigenvectors().transpose();
  Eigen::MatrixXd result = (Vm1 * Ssqrt).transpose();
  return result;
}

Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt(double etol) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());
  removedfunctions = 0;
  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  return es.eigenvectors() * diagonal.asDiagonal() *
         es.eigenvectors().transpose();
}

/***********************************
 * DIPOLE
 ***********************************/
void AODipole::Fill(const AOBasis& aobasis) {
  auto results = computeOneBodyIntegrals<libint2::Operator::emultipole1,
                                         std::array<libint2::Shell::real_t, 3>>(
      aobasis, _r);

  for (Index i = 0; i < 3; i++) {
    _aomatrix[i] = results[1 + i];  // emultipole1 returns: overlap, x-dipole,
                                    // y-dipole, z-dipole
  }
}
}  // namespace xtp
}  // namespace votca
