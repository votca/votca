/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/ERIs.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/threecenter.h"

// include libint last otherwise it overrides eigen
#include "votca/xtp/make_libint_work.h"
#define LIBINT2_CONSTEXPR_STATICS 0
#include <libint2.hpp>
#include <libint2/statics_definition.h>

namespace votca {
namespace xtp {

std::vector<std::vector<Index>> AOBasis::ComputeShellPairs(
    double threshold) const {

  Index nthreads = OPENMP::getMaxThreads();

  std::vector<libint2::Shell> shells = GenerateLibintBasis();

  // construct the 2-electron repulsion integrals engine
  std::vector<libint2::Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(libint2::Operator::overlap, getMaxNprim(), getMaxL(), 0);
  for (Index i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::vector<std::vector<Index>> pairs(shells.size());

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < Index(shells.size()); ++s1) {
    Index thread_id = OPENMP::getThreadId();

    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    Index n1 = shells[s1].size();

    for (Index s2 = 0; s2 <= s1; ++s2) {
      bool on_same_center = (shells[s1].O == shells[s2].O);
      bool significant = on_same_center;
      if (!on_same_center) {
        Index n2 = shells[s2].size();
        engine.compute(shells[s1], shells[s2]);
        Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
        significant = (buf_mat.norm() >= threshold);
      }
      if (significant) {
        pairs[s1].push_back(s2);
      }
    }
    std::sort(pairs[s1].begin(), pairs[s1].end());
  }
  return pairs;
}

using MatrixLibInt =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <libint2::Operator obtype,
          typename OperatorParams =
              typename libint2::operator_traits<obtype>::oper_params_type>
std::array<MatrixLibInt, libint2::operator_traits<obtype>::nopers>
    computeOneBodyIntegrals(const AOBasis& aobasis,
                            OperatorParams oparams = OperatorParams()) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();

  std::vector<std::vector<Index>> shellpair_list = aobasis.ComputeShellPairs();

  Index nopers = static_cast<Index>(libint2::operator_traits<obtype>::nopers);
  std::array<MatrixLibInt, libint2::operator_traits<obtype>::nopers> result;
  for (MatrixLibInt& r : result) {
    r = MatrixLibInt::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
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
        Eigen::Map<const MatrixLibInt> buf_mat(buf[op], n1, n2);
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
  aomatrix_ = computeOneBodyIntegrals<libint2::Operator::kinetic>(aobasis)[0];
}

/***********************************
 * OVERLAP
 ***********************************/
void AOOverlap::Fill(const AOBasis& aobasis) {
  aomatrix_ = computeOneBodyIntegrals<libint2::Operator::overlap>(aobasis)[0];
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

    Eigen::Map<const MatrixLibInt> buf_mat(buf[0], shell.getNumFunc(),
                                           shell.getNumFunc());
    return buf_mat;
  } catch (const libint2::Engine::lmax_exceeded& error) {
    std::ostringstream oss;
    oss << "\nA libint error occured:\n"
        << error.what() << "\n"
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

/***********************************
 * DIPOLE
 ***********************************/
void AODipole::Fill(const AOBasis& aobasis) {
  auto results = computeOneBodyIntegrals<libint2::Operator::emultipole1,
                                         std::array<libint2::Shell::real_t, 3>>(
      aobasis, r_);

  for (Index i = 0; i < 3; i++) {
    aomatrix_[i] = results[1 + i];  // emultipole1 returns: overlap, x-dipole,
                                    // y-dipole, z-dipole
  }
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

  aomatrix_ =
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

      Eigen::Map<const MatrixLibInt> buf_mat(buf[0], n1, n2);
      aomatrix_.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                       // {s2,s1} block, note the transpose!
        aomatrix_.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  }
}

Eigen::MatrixXd ERIs::ComputeSchwarzShells(const AOBasis& basis) const {

  Index noshells = basis.getNumofShells();

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(noshells, noshells);
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Engine> engines(nthreads);
  double epsilon = 0.0;
  engines[0] = libint2::Engine(libint2::Operator::coulomb, basis.getMaxNprim(),
                               static_cast<int>(basis.getMaxL()), 0, epsilon);

  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<libint2::Shell> shells = basis.GenerateLibintBasis();

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0l; s1 < basis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    Index n1 = shells[s1].size();

    for (Index s2 = 0l; s2 <= s1; ++s2) {
      Index n2 = shells[s2].size();
      Index n12 = n1 * n2;

      engines[thread_id]
          .compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              shells[s1], shells[s2], shells[s1], shells[s2]);

      Eigen::Map<const MatrixLibInt> buf_mat(buf[0], n12, n12);

      result(s2, s1) = std::sqrt(buf_mat.cwiseAbs().maxCoeff());
    }
  }
  return result.selfadjointView<Eigen::Upper>();
}

template <bool with_exchange>
std::array<Eigen::MatrixXd, 2> ERIs::Compute4c(const Eigen::MatrixXd& dmat,
                                               double error) const {
  assert(schwarzscreen_.rows() > 0 && schwarzscreen_.cols() > 0 &&
         "Please call Initialize_4c before running this");
  Index nthreads = OPENMP::getMaxThreads();

  Eigen::MatrixXd hartree = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  Eigen::MatrixXd exchange;
  if (with_exchange) {
    exchange = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  }
  Eigen::MatrixXd dnorm_block = ComputeShellBlockNorm(dmat);
  double fock_precision = error;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  Index max_nprim4 = maxnprim_ * maxnprim_ * maxnprim_ * maxnprim_;
  double engine_precision = std::min(fock_precision / dnorm_block.maxCoeff(),
                                     std::numeric_limits<double>::epsilon()) /
                            double(max_nprim4);
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(libint2::Operator::coulomb, int(maxnprim_),
                               int(maxL_), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }
  Index nshells = basis_.size();

#pragma omp parallel for schedule(dynamic)reduction(+ : hartree)reduction(+: exchange)
  for (Index s1 = 0; s1 < nshells; ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const auto& buf = engine.results();
    Index start_1 = starts_[s1];
    const libint2::Shell& shell1 = basis_[s1];
    Index n1 = shell1.size();

    auto sp12_iter = shellpairdata_[s1].begin();
    for (Index s2 : shellpairs_[s1]) {
      Index start_2 = starts_[s2];
      const libint2::Shell& shell2 = basis_[s2];
      Index n2 = shell2.size();
      double dnorm_12 = dnorm_block(s1, s2);
      const libint2::ShellPair* sp12 = &(*sp12_iter);
      ++sp12_iter;

      for (Index s3 = 0; s3 <= s1; ++s3) {

        Index start_3 = starts_[s3];
        const libint2::Shell& shell3 = basis_[s3];
        Index n3 = shell3.size();
        auto sp34_iter = shellpairdata_[s3].begin();
        double dnorm_123 = std::max(dnorm_block(s1, s3),
                                    std::max(dnorm_block(s2, s3), dnorm_12));
        Index s4max = (s1 == s3) ? s2 : s3;
        for (Index s4 : shellpairs_[s3]) {
          if (s4 > s4max) {
            break;
          }  // for each s3, s4 are stored in monotonically increasing
             // order

          const libint2::ShellPair* sp34 = &(*sp34_iter);
          // must update the iter even if going to skip s4
          ++sp34_iter;
          double dnorm_1234 =
              std::max(dnorm_block(s1, s4),
                       std::max(dnorm_block(s2, s4),
                                std::max(dnorm_block(s3, s4), dnorm_123)));

          if (dnorm_1234 * schwarzscreen_(s1, s2) * schwarzscreen_(s3, s4) <
              fock_precision) {
            continue;
          }

          const libint2::Shell& shell4 = basis_[s4];
          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                  shell1, shell2, shell3, shell4, sp12, sp34);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr) {
            continue;  // if all integrals screened out, skip to next quartet
          }
          Index start_4 = starts_[s4];
          Index n4 = shell4.size();
          Index s12_deg = (s1 == s2) ? 1 : 2;
          Index s34_deg = (s3 == s4) ? 1 : 2;
          Index s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1 : 2) : 2;
          Index s1234_deg = s12_deg * s34_deg * s12_34_deg;

          for (Index f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const Index bf1 = f1 + start_1;
            for (Index f2 = 0; f2 != n2; ++f2) {
              const Index bf2 = f2 + start_2;
              for (Index f3 = 0; f3 != n3; ++f3) {
                const Index bf3 = f3 + start_3;
                for (Index f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const Index bf4 = f4 + start_4;

                  const double value = buf_1234[f1234];

                  const double value_scal_by_deg = value * double(s1234_deg);

                  hartree(bf1, bf2) += dmat(bf3, bf4) * value_scal_by_deg;
                  hartree(bf3, bf4) += dmat(bf1, bf2) * value_scal_by_deg;
                  if (with_exchange) {
                    exchange(bf1, bf3) -= dmat(bf2, bf4) * value_scal_by_deg;
                    exchange(bf2, bf3) -= dmat(bf1, bf4) * value_scal_by_deg;
                    exchange(bf2, bf4) -= dmat(bf1, bf3) * value_scal_by_deg;
                    exchange(bf1, bf4) -= dmat(bf2, bf3) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::array<Eigen::MatrixXd, 2> result2;
  // 0.25=0.5(symmetrisation)*0.5(our dmat has a factor 2)
  result2[0] = 0.25 * (hartree + hartree.transpose());
  if (with_exchange) {
    // prefactor
    result2[1] = 0.125 * (exchange + exchange.transpose());
  }
  return result2;
}

template std::array<Eigen::MatrixXd, 2> ERIs::Compute4c<true>(
    const Eigen::MatrixXd& dmat, double error) const;
template std::array<Eigen::MatrixXd, 2> ERIs::Compute4c<false>(
    const Eigen::MatrixXd& dmat, double error) const;

void TCMatrix_dft::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis) {
  {
    AOCoulomb auxAOcoulomb;
    auxAOcoulomb.Fill(auxbasis);
    inv_sqrt_ = auxAOcoulomb.Pseudo_InvSqrt(1e-8);
    removedfunctions_ = auxAOcoulomb.Removedfunctions();
  }
  matrix_ = std::vector<Symmetric_Matrix>(auxbasis.AOBasisSize());

#pragma omp parallel for schedule(dynamic, 4)
  for (Index i = 0; i < auxbasis.AOBasisSize(); i++) {
    matrix_[i] = Symmetric_Matrix(dftbasis.AOBasisSize());
  }

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

#pragma omp parallel for schedule(dynamic)
  for (Index is = dftbasis.getNumofShells() - 1; is >= 0; is--) {

    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    const libint2::Shell& dftshell = dftshells[is];
    Index start = shell2bf[is];
    std::vector<Eigen::MatrixXd> block(dftshell.size());
    for (Index i = 0; i < Index(dftshell.size()); i++) {
      Index size = start + i + 1;
      block[i] = Eigen::MatrixXd::Zero(auxbasis.AOBasisSize(), size);
    }

    for (Index aux = 0; aux < auxbasis.getNumofShells(); aux++) {
      const libint2::Shell& auxshell = auxshells[aux];
      Index aux_start = auxshell2bf[aux];

      for (Index dis = 0; dis <= is; dis++) {

        const libint2::Shell& shell_col = dftshells[dis];
        Index col_start = shell2bf[dis];
        engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
            auxshell, libint2::Shell::unit(), dftshell, shell_col);

        if (buf[0] == nullptr) {
          continue;
        }
        Eigen::TensorMap<Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
            result(buf[0], auxshell.size(), dftshell.size(), shell_col.size());

        for (size_t left = 0; left < dftshell.size(); left++) {
          for (size_t auxf = 0; auxf < auxshell.size(); auxf++) {
            for (size_t col = 0; col < shell_col.size(); col++) {
              // symmetry
              if ((col_start + col) > (start + left)) {
                break;
              }
              block[left](aux_start + auxf, col_start + col) =
                  result(auxf, left, col);
            }
          }
        }
      }
    }

    for (Index i = 0; i < Index(block.size()); ++i) {
      Eigen::MatrixXd temp = inv_sqrt_ * block[i];
      for (Index mu = 0; mu < temp.rows(); ++mu) {
        for (Index j = 0; j < temp.cols(); ++j) {
          matrix_[mu](i + start, j) = temp(mu, j);
        }
      }
    }
  }

  return;
}

/*
 * Determines the 3-center integrals for a given shell in the aux basis
 * by calculating the 3-center repulsion integral of the functions in the
 * aux shell with ALL functions in the DFT basis set
 */
std::vector<Eigen::MatrixXd> ComputeAO3cBlock(const libint2::Shell& auxshell,
                                              const AOBasis& dftbasis,
                                              libint2::Engine& engine) {
  std::vector<Eigen::MatrixXd> ao3c = std::vector<Eigen::MatrixXd>(
      auxshell.size(),
      Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));

  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();

  const libint2::Engine::target_ptr_vec& buf = engine.results();
  // alpha-loop over the "left" DFT basis function
  for (Index row = 0; row < Index(dftshells.size()); row++) {

    const libint2::Shell& shell_row = dftshells[row];
    const Index row_start = shell2bf[row];
    // ThreecMatrix is symmetric, restrict explicit calculation to triangular
    // matrix
    for (Index col = 0; col <= row; col++) {
      const libint2::Shell& shell_col = dftshells[col];
      const Index col_start = shell2bf[col];

      engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
          auxshell, libint2::Shell::unit(), shell_col, shell_row);

      if (buf[0] == nullptr) {
        continue;
      }
      Eigen::TensorMap<Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
          result(buf[0], auxshell.size(), shell_col.size(), shell_row.size());

      for (size_t aux_c = 0; aux_c < auxshell.size(); aux_c++) {
        for (size_t col_c = 0; col_c < shell_col.size(); col_c++) {
          for (size_t row_c = 0; row_c < shell_row.size(); row_c++) {
            // ao3c is col major result is row-major so this is ideal
            ao3c[aux_c](row_start + row_c, col_start + col_c) =
                result(aux_c, col_c, row_c);
          }  // ROW copy
        }    // COL copy
      }      // AUX copy

    }  // gamma-loop
  }    // alpha-loop

  for (Eigen::MatrixXd& mat : ao3c) {
    mat.triangularView<Eigen::Upper>() =
        mat.triangularView<Eigen::Lower>().transpose();
  }
  return ao3c;
}

void TCMatrix_gwbse::Fill3cMO(const AOBasis& auxbasis, const AOBasis& dftbasis,
                              const Eigen::MatrixXd& dft_orbitals) {

  const Eigen::MatrixXd dftm = dft_orbitals.middleCols(mmin_, mtotal_);
  const Eigen::MatrixXd dftn =
      dft_orbitals.middleCols(nmin_, ntotal_).transpose();

  OpenMP_CUDA transform;
  transform.setOperators(dftn, dftm);
  Index nthreads = OPENMP::getMaxThreads();

  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index aux = 0; aux < Index(auxshells.size()); aux++) {
      const libint2::Shell& auxshell = auxshells[aux];

      std::vector<Eigen::MatrixXd> ao3c =
          ComputeAO3cBlock(auxshell, dftbasis, engines[threadid]);

      // this is basically a transpose of AO3c and at the same time the ao->mo
      // transformation
      // we do not want to put it into  matrix_ straight away is because,
      //  matrix_ is shared between all threads and we want a nice clean access
      // pattern to it
      std::vector<Eigen::MatrixXd> block = std::vector<Eigen::MatrixXd>(
          mtotal_, Eigen::MatrixXd::Zero(ntotal_, ao3c.size()));

      Index dim = static_cast<Index>(ao3c.size());
      for (Index k = 0; k < dim; ++k) {
        transform.MultiplyLeftRight(ao3c[k], threadid);
        for (Index i = 0; i < ao3c[k].cols(); ++i) {
          block[i].col(k) = ao3c[k].col(i);
        }
      }

      // put into correct position
      for (Index m_level = 0; m_level < mtotal_; m_level++) {
        matrix_[m_level].middleCols(auxshell2bf[aux], auxshell.size()) =
            block[m_level];
      }  // m-th DFT orbital
    }    // shells of GW basis set
  }
}

}  // namespace xtp
}  // namespace votca
