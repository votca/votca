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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/ERIs.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/symmetric_matrix.h>

namespace votca {
namespace xtp {

void ERIs::Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis) {
  _threecenter.Fill(auxbasis, dftbasis);
  return;
}

void ERIs::Initialize_4c_small_molecule(const AOBasis& dftbasis) {
  _fourcenter.Fill_4c_small_molecule(dftbasis);
  return;
}

void ERIs::Initialize_4c_screening(const AOBasis& dftbasis, double eps) {
  _with_screening = true;
  _screening_eps = eps;
  CalculateERIsDiagonals(dftbasis);
  return;
}

Mat_p_Energy ERIs::CalculateERIs(const Eigen::MatrixXd& DMAT) const {
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> ERIS_thread = std::vector<Eigen::MatrixXd>(
      nthreads, Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()));
  Symmetric_Matrix dmat_sym = Symmetric_Matrix(DMAT);
#pragma omp parallel for
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Symmetric_Matrix& threecenter = _threecenter[i];
    // Trace over prod::DMAT,I(l)=componentwise product over
    const double factor = threecenter.TraceofProd(dmat_sym);
    Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper> m =
        ERIS_thread[OPENMP::getThreadId()].selfadjointView<Eigen::Upper>();
    threecenter.AddtoEigenUpperMatrix(m, factor);
  }

  Eigen::MatrixXd ERIs2 =
      std::accumulate(ERIS_thread.begin(), ERIS_thread.end(),
                      Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()).eval());
  ERIs2 = ERIs2.selfadjointView<Eigen::Upper>();
  double energy = CalculateEnergy(DMAT, ERIs2);
  return Mat_p_Energy(energy, ERIs2);
}

Eigen::MatrixXcd ERIs::ContractRightIndecesWithMatrix(const Eigen::MatrixXcd& mat) const {
  int nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXcd> ERIS_thread = std::vector<Eigen::MatrixXcd>(
      nthreads, Eigen::MatrixXcd::Zero(mat.rows(), mat.cols()));
#pragma omp parallel for
  for (int i = 0; i < _threecenter.size(); i++) {
    const Symmetric_Matrix& threecenter = _threecenter[i];
    // Trace over prod::DMAT,I(l)=componentwise product over
    const std::complex<double> factor = (threecenter.FullMatrix().cwiseProduct(mat)).sum();
    threecenter.AddtoEigenMatrix(ERIS_thread[OPENMP::getThreadId()], factor);
  }

  Eigen::MatrixXcd ERIs =
      std::accumulate(ERIS_thread.begin(), ERIS_thread.end(),
                      Eigen::MatrixXcd::Zero(mat.rows(), mat.cols()).eval());
  return ERIs;
}


Mat_p_Energy ERIs::CalculateEXX(const Eigen::MatrixXd& DMAT) const {

  std::vector<Eigen::MatrixXd> EXX_thread = std::vector<Eigen::MatrixXd>(
      OPENMP::getMaxThreads(), Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()));

#pragma omp parallel for
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Eigen::MatrixXd threecenter = _threecenter[i].UpperMatrix();
    EXX_thread[OPENMP::getThreadId()] +=
        threecenter.selfadjointView<Eigen::Upper>() * DMAT *
        threecenter.selfadjointView<Eigen::Upper>();
  }
  Eigen::MatrixXd EXXs =
      std::accumulate(EXX_thread.begin(), EXX_thread.end(),
                      Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()).eval());
  double energy = CalculateEnergy(DMAT, EXXs);
  return Mat_p_Energy(energy, EXXs);
}

Mat_p_Energy ERIs::CalculateEXX(const Eigen::MatrixXd& occMos,
                                const Eigen::MatrixXd& DMAT) const {

  std::vector<Eigen::MatrixXd> EXX_thread = std::vector<Eigen::MatrixXd>(
      OPENMP::getMaxThreads(),
      Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows()));

#pragma omp parallel for
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Eigen::MatrixXd TCxMOs_T =
        occMos.transpose() *
        _threecenter[i].UpperMatrix().selfadjointView<Eigen::Upper>();
    EXX_thread[OPENMP::getThreadId()] += TCxMOs_T.transpose() * TCxMOs_T;
  }

  Eigen::MatrixXd EXXs =
      2 * std::accumulate(
              EXX_thread.begin(), EXX_thread.end(),
              Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows()).eval());

  double energy = CalculateEnergy(DMAT, EXXs);
  return Mat_p_Energy(energy, EXXs);
}

Mat_p_Energy ERIs::CalculateERIs_4c_small_molecule(
    const Eigen::MatrixXd& DMAT) const {

  Eigen::MatrixXd ERIs2 = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

  const Eigen::VectorXd& fourc_vector = _fourcenter.get_4c_vector();

  Index dftBasisSize = DMAT.rows();
  Index vectorSize = (dftBasisSize * (dftBasisSize + 1)) / 2;
#pragma omp parallel for
  for (Index i = 0; i < dftBasisSize; i++) {
    Index sum_i = (i * (i + 1)) / 2;
    for (Index j = i; j < dftBasisSize; j++) {
      Index index_ij = dftBasisSize * i - sum_i + j;
      Index index_ij_kl_a =
          vectorSize * index_ij - (index_ij * (index_ij + 1)) / 2;
      for (Index k = 0; k < dftBasisSize; k++) {
        Index sum_k = (k * (k + 1)) / 2;
        for (Index l = k; l < dftBasisSize; l++) {
          Index index_kl = dftBasisSize * k - sum_k + l;

          Index index_ij_kl = index_ij_kl_a + index_kl;
          if (index_ij > index_kl) {
            index_ij_kl = vectorSize * index_kl -
                          (index_kl * (index_kl + 1)) / 2 + index_ij;
          }

          if (l == k) {
            ERIs2(i, j) += DMAT(k, l) * fourc_vector(index_ij_kl);
          } else {
            ERIs2(i, j) += 2. * DMAT(k, l) * fourc_vector(index_ij_kl);
          }
        }
      }
      ERIs2(j, i) = ERIs2(i, j);
    }
  }

  double energy = CalculateEnergy(DMAT, ERIs2);
  return Mat_p_Energy(energy, ERIs2);
}

Mat_p_Energy ERIs::CalculateEXX_4c_small_molecule(
    const Eigen::MatrixXd& DMAT) const {
  std::vector<Eigen::MatrixXd> EXX_thread = std::vector<Eigen::MatrixXd>(
      OPENMP::getMaxThreads(), Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()));

  const Eigen::VectorXd& fourc_vector = _fourcenter.get_4c_vector();

  Index dftBasisSize = DMAT.rows();
  Index vectorSize = (dftBasisSize * (dftBasisSize + 1)) / 2;
#pragma omp parallel for
  for (Index i = 0; i < dftBasisSize; i++) {
    Index thread = OPENMP::getThreadId();
    Index sum_i = (i * (i + 1)) / 2;
    for (Index j = i; j < dftBasisSize; j++) {
      Index index_ij = DMAT.cols() * i - sum_i + j;
      Index index_ij_kl_a =
          vectorSize * index_ij - (index_ij * (index_ij + 1)) / 2;
      for (Index k = 0; k < dftBasisSize; k++) {
        Index sum_k = (k * (k + 1)) / 2;
        for (Index l = k; l < dftBasisSize; l++) {
          Index index_kl = DMAT.cols() * k - sum_k + l;

          Index _index_ij_kl = index_ij_kl_a + index_kl;
          if (index_ij > index_kl) {
            _index_ij_kl = vectorSize * index_kl -
                           (index_kl * (index_kl + 1)) / 2 + index_ij;
          }
          double factorij = 1;
          if (i == j) {
            factorij = 0.5;
          }
          double factorkl = 1;
          if (l == k) {
            factorkl = 0.5;
          }
          double factor = factorij * factorkl;
          EXX_thread[thread](i, l) +=
              factor * DMAT(j, k) * fourc_vector(_index_ij_kl);
          EXX_thread[thread](j, l) +=
              factor * DMAT(i, k) * fourc_vector(_index_ij_kl);
          EXX_thread[thread](i, k) +=
              factor * DMAT(j, l) * fourc_vector(_index_ij_kl);
          EXX_thread[thread](j, k) +=
              factor * DMAT(i, l) * fourc_vector(_index_ij_kl);
        }
      }
    }
  }

  Eigen::MatrixXd EXXs =
      std::accumulate(EXX_thread.begin(), EXX_thread.end(),
                      Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols()).eval());

  double energy = CalculateEnergy(DMAT, EXXs);
  return Mat_p_Energy(energy, EXXs);
}

Mat_p_Energy ERIs::CalculateERIs_4c_direct(const AOBasis& dftbasis,
                                           const Eigen::MatrixXd& DMAT) const {

  // Number of shells
  Index numShells = dftbasis.getNumofShells();

  // Initialize ERIs matrix
  Eigen::MatrixXd ERIs2 = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

#pragma omp parallel
  {  // Begin omp parallel

    Eigen::MatrixXd ERIs_thread =
        Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

#pragma omp for
    for (Index iShell_3 = 0; iShell_3 < numShells; iShell_3++) {
      const AOShell& shell_3 = dftbasis.getShell(iShell_3);
      Index numFunc_3 = shell_3.getNumFunc();
      for (Index iShell_4 = iShell_3; iShell_4 < numShells; iShell_4++) {
        const AOShell& shell_4 = dftbasis.getShell(iShell_4);
        Index numFunc_4 = shell_4.getNumFunc();
        for (Index iShell_1 = iShell_3; iShell_1 < numShells; iShell_1++) {
          const AOShell& shell_1 = dftbasis.getShell(iShell_1);
          Index numFunc_1 = shell_1.getNumFunc();
          for (Index iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
            const AOShell& shell_2 = dftbasis.getShell(iShell_2);
            Index numFunc_2 = shell_2.getNumFunc();

            // Pre-screening
            if (_with_screening && CheckScreen(_screening_eps, shell_1, shell_2,
                                               shell_3, shell_4)) {
              continue;
            }

            // Get the current 4c block
            Eigen::Tensor<double, 4> block(numFunc_1, numFunc_2, numFunc_3,
                                           numFunc_4);
            block.setZero();
            bool nonzero = _fourcenter.FillFourCenterRepBlock(
                block, shell_1, shell_2, shell_3, shell_4);

            // If there are only zeros, we don't need to put anything in the
            // ERIs matrix
            if (!nonzero) {
              continue;
            }

            // Begin fill ERIs matrix

            FillERIsBlock<false>(ERIs_thread, DMAT, block, shell_1, shell_2,
                                 shell_3, shell_4);

            // Symmetry 1 <--> 2
            if (iShell_1 != iShell_2) {
              FillERIsBlock<false>(ERIs_thread, DMAT, block, shell_2, shell_1,
                                   shell_3, shell_4);
            }

            // Symmetry 3 <--> 4
            if (iShell_3 != iShell_4) {
              FillERIsBlock<false>(ERIs_thread, DMAT, block, shell_1, shell_2,
                                   shell_4, shell_3);
            }

            // Symmetry 1 <--> 2 and 3 <--> 4
            if (iShell_1 != iShell_2 && iShell_3 != iShell_4) {
              FillERIsBlock<false>(ERIs_thread, DMAT, block, shell_2, shell_1,
                                   shell_4, shell_3);
            }

            // Symmetry (1, 2) <--> (3, 4)
            if (iShell_1 != iShell_3) {

              FillERIsBlock<true>(ERIs_thread, DMAT, block, shell_3, shell_4,
                                  shell_1, shell_2);

              // Symmetry 1 <--> 2
              if (iShell_1 != iShell_2) {
                FillERIsBlock<true>(ERIs_thread, DMAT, block, shell_3, shell_4,
                                    shell_2, shell_1);
              }

              // Symmetry 3 <--> 4
              if (iShell_3 != iShell_4) {
                FillERIsBlock<true>(ERIs_thread, DMAT, block, shell_4, shell_3,
                                    shell_1, shell_2);
              }

              // Symmetry 1 <--> 2 and 3 <--> 4
              if (iShell_1 != iShell_2 && iShell_3 != iShell_4) {
                FillERIsBlock<true>(ERIs_thread, DMAT, block, shell_4, shell_3,
                                    shell_2, shell_1);
              }
            }

            // End fill ERIs matrix
          }  // End loop over shell 2
        }    // End loop over shell 1
      }      // End loop over shell 4
    }        // End loop over shell 3

#pragma omp critical
    { ERIs2 += ERIs_thread; }
  }

  ERIs2 = ERIs2.selfadjointView<Eigen::Upper>();
  double energy = CalculateEnergy(DMAT, ERIs2);
  return Mat_p_Energy(energy, ERIs2);
}

template <bool transposed_block>
void ERIs::FillERIsBlock(Eigen::MatrixXd& ERIsCur, const Eigen::MatrixXd& DMAT,
                         const Eigen::Tensor<double, 4>& block,
                         const AOShell& shell_1, const AOShell& shell_2,
                         const AOShell& shell_3, const AOShell& shell_4) const {

  for (Index iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
    Index ind_3 = shell_3.getStartIndex() + iFunc_3;
    for (Index iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
      Index ind_4 = shell_4.getStartIndex() + iFunc_4;

      // Symmetry
      if (ind_3 > ind_4) {
        continue;
      }

      for (Index iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
        Index ind_1 = shell_1.getStartIndex() + iFunc_1;
        for (Index iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
          Index ind_2 = shell_2.getStartIndex() + iFunc_2;

          // Symmetry
          if (ind_1 > ind_2) {
            continue;
          }

          // Symmetry for diagonal elements
          double multiplier = (ind_1 == ind_2 ? 1.0 : 2.0);
          // Fill ERIs matrix
          if (!transposed_block) {
            ERIsCur(ind_3, ind_4) += multiplier * DMAT(ind_1, ind_2) *
                                     block(iFunc_1, iFunc_2, iFunc_3, iFunc_4);
          } else {
            ERIsCur(ind_3, ind_4) += multiplier * DMAT(ind_1, ind_2) *
                                     block(iFunc_3, iFunc_4, iFunc_1, iFunc_2);
          }

        }  // End loop over functions in shell 2
      }    // End loop over functions in shell 1
    }      // End loop over functions in shell 4
  }        // End loop over functions in shell 3

  return;
}

void ERIs::CalculateERIsDiagonals(const AOBasis& dftbasis) {
  // Number of shells
  Index numShells = dftbasis.getNumofShells();
  // Total number of functions
  Index dftBasisSize = dftbasis.AOBasisSize();

  _diagonals = Eigen::MatrixXd::Zero(dftBasisSize, dftBasisSize);

  for (Index iShell_1 = 0; iShell_1 < numShells; iShell_1++) {
    const AOShell& shell_1 = dftbasis.getShell(iShell_1);
    Index numFunc_1 = shell_1.getNumFunc();
    for (Index iShell_2 = iShell_1; iShell_2 < numShells; iShell_2++) {
      const AOShell& shell_2 = dftbasis.getShell(iShell_2);
      Index numFunc_2 = shell_2.getNumFunc();

      // Get the current 4c block
      Eigen::Tensor<double, 4> block(numFunc_1, numFunc_2, numFunc_1,
                                     numFunc_2);
      block.setZero();
      bool nonzero = _fourcenter.FillFourCenterRepBlock(block, shell_1, shell_2,
                                                        shell_1, shell_2);

      if (!nonzero) {
        continue;
      }

      for (Index iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
        Index ind_1 = shell_1.getStartIndex() + iFunc_1;
        for (Index iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
          Index ind_2 = shell_2.getStartIndex() + iFunc_2;

          // Symmetry
          if (ind_1 > ind_2) {
            continue;
          }

          _diagonals(ind_1, ind_2) = block(iFunc_1, iFunc_2, iFunc_1, iFunc_2);

          // Symmetry
          if (ind_1 != ind_2) {
            _diagonals(ind_2, ind_1) = _diagonals(ind_1, ind_2);
          }
        }
      }
    }
  }

  return;
}

bool ERIs::CheckScreen(double eps, const AOShell& shell_1,
                       const AOShell& shell_2, const AOShell& shell_3,
                       const AOShell& shell_4) const {

  const double eps2 = eps * eps;

  for (Index iFunc_3 = 0; iFunc_3 < shell_3.getNumFunc(); iFunc_3++) {
    Index ind_3 = shell_3.getStartIndex() + iFunc_3;
    for (Index iFunc_4 = 0; iFunc_4 < shell_4.getNumFunc(); iFunc_4++) {
      Index ind_4 = shell_4.getStartIndex() + iFunc_4;

      // Symmetry
      if (ind_3 > ind_4) {
        continue;
      }

      for (Index iFunc_1 = 0; iFunc_1 < shell_1.getNumFunc(); iFunc_1++) {
        Index ind_1 = shell_1.getStartIndex() + iFunc_1;
        for (Index iFunc_2 = 0; iFunc_2 < shell_2.getNumFunc(); iFunc_2++) {
          Index ind_2 = shell_2.getStartIndex() + iFunc_2;

          // Symmetry
          if (ind_1 > ind_2) {
            continue;
          }

          // Cauchy–Schwarz
          // <ab|cd> <= sqrt(<ab|ab>) * sqrt(<cd|cd>)
          double ub = _diagonals(ind_1, ind_2) * _diagonals(ind_3, ind_4);

          // Compare with tolerance
          if (ub > eps2) {
            return false;  // We must compute ERIS for the whole block
          }
        }
      }
    }
  }

  return true;  // We can skip the whole block
}

double ERIs::CalculateEnergy(const Eigen::MatrixXd& DMAT,
                             const Eigen::MatrixXd& matrix_operator) const {
  return matrix_operator.cwiseProduct(DMAT).sum();
}

}  // namespace xtp
}  // namespace votca
