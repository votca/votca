/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE dftengine_private_test
#include "xtp_libint2.h"
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <votca/xtp/basisset.h>

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/dftengine.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

namespace votca {
namespace xtp {

constexpr double kTol = 1e-10;

class DFTEngineTestAccess {
 public:
  static double NuclearRepulsion(const DFTEngine& e, const QMMolecule& mol) {
    return e.NuclearRepulsion(mol);
  }

  static Eigen::MatrixXd SphericalAverageShells(const DFTEngine& e,
                                                const Eigen::MatrixXd& dmat,
                                                const AOBasis& basis) {
    return e.SphericalAverageShells(dmat, basis);
  }

  static Eigen::VectorXd BuildEHTOrbitalEnergies(const DFTEngine& e,
                                                 const QMMolecule& mol) {
    return e.BuildEHTOrbitalEnergies(mol);
  }

  static Eigen::MatrixXd BuildEHTHamiltonian(const DFTEngine& e,
                                             const QMMolecule& mol) {
    return e.BuildEHTHamiltonian(mol);
  }

  static Eigen::MatrixXd OrthogonalizeGuess(const DFTEngine& e,
                                            const Eigen::MatrixXd& guess) {
    return e.OrthogonalizeGuess(guess);
  }

  static Eigen::MatrixXd InsertZeroRows(DFTEngine& e, Eigen::MatrixXd M,
                                        Index startidx, Index numofzerorows) {
    return e.InsertZeroRows(M, startidx, numofzerorows);
  }

  static Eigen::MatrixXd InsertZeroCols(DFTEngine& e, Eigen::MatrixXd M,
                                        Index startidx, Index numofzerocols) {
    return e.InsertZeroCols(M, startidx, numofzerocols);
  }

  static void SetBasis(DFTEngine& e, const AOBasis& basis) {
    e.dftbasis_ = basis;
  }

  static void FillOverlap(DFTEngine& e) { e.dftAOoverlap_.Fill(e.dftbasis_); }

  static void SetBasisName(DFTEngine& e, const std::string& basis_name) {
    e.dftbasis_name_ = basis_name;
  }

  static const AOOverlap& Overlap(const DFTEngine& e) {
    return e.dftAOoverlap_;
  }
};

QMMolecule MakeH2(double distance_bohr) {
  QMMolecule mol("H2", 0);
  mol.push_back(QMAtom(0, "H", Eigen::Vector3d(0.0, 0.0, 0.0)));
  mol.push_back(QMAtom(1, "H", Eigen::Vector3d(distance_bohr, 0.0, 0.0)));
  return mol;
}

QMMolecule MakeSingleAtom(const std::string& element) {
  QMMolecule mol("single_atom", 0);
  mol.push_back(QMAtom(0, element, Eigen::Vector3d::Zero()));
  return mol;
}

AOBasis MakeBasis(const std::string& basis_name, const QMMolecule& mol) {
  BasisSet basisset;
  basisset.Load(basis_name);

  AOBasis basis;
  basis.Fill(basisset, mol);
  return basis;
}

bool BlockIsConstant(const Eigen::MatrixXd& m, Index row0, Index col0,
                     Index rows, Index cols, double value, double tol = 1e-12) {
  return ((m.block(row0, col0, rows, cols).array() - value).abs() < tol).all();
}

}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_CASE(nuclear_repulsion_h2_at_2_bohr) {
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  QMMolecule mol = MakeH2(2.0);

  // Z1 Z2 / R = 1 * 1 / 2 = 0.5 Ha
  const double e = DFTEngineTestAccess::NuclearRepulsion(engine, mol);

  BOOST_TEST(e == 0.5, boost::test_tools::tolerance(kTol));
}

BOOST_AUTO_TEST_CASE(build_eht_hamiltonian_is_symmetric_and_matches_formula) {
  libint2::initialize();
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  const std::string basis_name = "3-21G";
  QMMolecule mol = MakeH2(1.4);

  AOBasis basis = MakeBasis(basis_name, mol);
  DFTEngineTestAccess::SetBasis(engine, basis);
  DFTEngineTestAccess::SetBasisName(engine, basis_name);
  DFTEngineTestAccess::FillOverlap(engine);

  const Eigen::MatrixXd H =
      DFTEngineTestAccess::BuildEHTHamiltonian(engine, mol);
  const Eigen::VectorXd eps =
      DFTEngineTestAccess::BuildEHTOrbitalEnergies(engine, mol);
  const Eigen::MatrixXd& S = DFTEngineTestAccess::Overlap(engine).Matrix();

  BOOST_REQUIRE_EQUAL(H.rows(), basis.AOBasisSize());
  BOOST_REQUIRE_EQUAL(H.cols(), basis.AOBasisSize());

  BOOST_TEST((H - H.transpose()).norm() < 1e-12);

  for (Index mu = 0; mu < basis.AOBasisSize(); ++mu) {
    BOOST_TEST(H(mu, mu) == eps(mu), boost::test_tools::tolerance(kTol));
  }

  constexpr double K = 1.75;
  for (Index mu = 0; mu < basis.AOBasisSize(); ++mu) {
    for (Index nu = 0; nu < mu; ++nu) {
      const double expected = K * S(mu, nu) * 0.5 * (eps(mu) + eps(nu));
      BOOST_TEST(H(mu, nu) == expected, boost::test_tools::tolerance(1e-12));
      BOOST_TEST(H(nu, mu) == expected, boost::test_tools::tolerance(1e-12));
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(spherical_average_shells_makes_shell_blocks_uniform) {
  libint2::initialize();
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  // Needs at least one multi-function shell; O/STO-3G is a reasonable choice.
  const std::string basis_name = "3-21G";
  QMMolecule mol = MakeSingleAtom("O");
  AOBasis basis = MakeBasis(basis_name, mol);

  const Index nao = basis.AOBasisSize();
  BOOST_REQUIRE(nao > 1);

  Eigen::MatrixXd dmat = Eigen::MatrixXd::Zero(nao, nao);
  for (Index i = 0; i < nao; ++i) {
    for (Index j = 0; j < nao; ++j) {
      dmat(i, j) =
          10.0 * static_cast<double>(i + 1) + static_cast<double>(j + 1);
    }
  }

  const Eigen::MatrixXd avg =
      DFTEngineTestAccess::SphericalAverageShells(engine, dmat, basis);

  BOOST_REQUIRE_EQUAL(avg.rows(), nao);
  BOOST_REQUIRE_EQUAL(avg.cols(), nao);

  for (const AOShell& shellrow : basis) {
    const Index size_row = shellrow.getNumFunc();
    const Index start_row = shellrow.getStartIndex();

    for (const AOShell& shellcol : basis) {
      const Index size_col = shellcol.getNumFunc();
      const Index start_col = shellcol.getStartIndex();

      const Eigen::MatrixXd in_block =
          dmat.block(start_row, start_col, size_row, size_col);
      const Eigen::MatrixXd out_block =
          avg.block(start_row, start_col, size_row, size_col);

      if (shellrow.getL() == shellcol.getL()) {
        const double diagavg =
            in_block.diagonal().sum() / static_cast<double>(in_block.rows());

        if (size_row == 1 && size_col == 1) {
          BOOST_TEST(out_block(0, 0) == diagavg,
                     boost::test_tools::tolerance(1e-12));
        } else {
          const Index offdiag_n =
              in_block.rows() * in_block.cols() - in_block.cols();
          const double offdiagavg =
              (in_block.sum() - in_block.diagonal().sum()) /
              static_cast<double>(offdiag_n);

          for (Index i = 0; i < size_row; ++i) {
            for (Index j = 0; j < size_col; ++j) {
              if (i == j) {
                BOOST_TEST(out_block(i, j) == diagavg,
                           boost::test_tools::tolerance(1e-12));
              } else {
                BOOST_TEST(out_block(i, j) == offdiagavg,
                           boost::test_tools::tolerance(1e-12));
              }
            }
          }
        }
      } else {
        const double expected =
            in_block.sum() / static_cast<double>(in_block.size());
        BOOST_TEST(BlockIsConstant(avg, start_row, start_col, size_row,
                                   size_col, expected));
      }
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(orthogonalize_guess_produces_s_orthonormal_vectors) {
  libint2::initialize();

  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  const std::string basis_name = "3-21G";
  QMMolecule mol = MakeH2(1.4);

  AOBasis basis = MakeBasis(basis_name, mol);
  DFTEngineTestAccess::SetBasis(engine, basis);
  DFTEngineTestAccess::SetBasisName(engine, basis_name);
  DFTEngineTestAccess::FillOverlap(engine);

  const Index nao = basis.AOBasisSize();
  BOOST_REQUIRE(nao >= 2);

  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(nao, 2);
  guess(0, 0) = 1.0;
  guess(0, 1) = 1.0;
  guess(1, 0) = 0.25;
  guess(1, 1) = 1.0;
  if (nao > 2) {
    guess(2, 0) = -0.3;
    guess(2, 1) = 0.1;
  }

  const Eigen::MatrixXd orth =
      DFTEngineTestAccess::OrthogonalizeGuess(engine, guess);
  const Eigen::MatrixXd& S = DFTEngineTestAccess::Overlap(engine).Matrix();

  const Eigen::MatrixXd metric = orth.transpose() * S * orth;
  const Eigen::MatrixXd I =
      Eigen::MatrixXd::Identity(metric.rows(), metric.cols());

  BOOST_TEST((metric - I).norm() < 1e-10);

  libint2::finalize();
}
