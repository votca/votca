/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#define BOOST_TEST_MODULE gridboxaccumulator_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/gridboxaccumulator.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/vxc_grid.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(gridboxaccumulator_test)

namespace {
// a deterministic box block, different for every box and matrix
Eigen::MatrixXd Block(Index box, Index k, Index n) {
  Eigen::MatrixXd m(n, n);
  for (Index i = 0; i < n; ++i) {
    for (Index j = 0; j < n; ++j) {
      m(i, j) = std::sin(double(1 + box) * 0.37 + double(k) * 1.3 +
                         double(i) * 0.11 - double(j) * 0.07);
    }
  }
  return m;
}

// sums the blocks of all boxes into n matrices, with the accumulator in
// parallel (budget selects thread-local copies or per-matrix locks) or, for
// budget < 0, serially without it
std::vector<Eigen::MatrixXd> Accumulate(const Vxc_Grid& grid, Index n_ao,
                                        Index n, double budget) {
  std::vector<Eigen::MatrixXd> result(std::size_t(n),
                                      Eigen::MatrixXd::Zero(n_ao, n_ao));
  if (budget < 0) {
    for (Index i = 0; i < grid.getBoxesSize(); ++i) {
      const GridBox& box = grid[i];
      if (!box.Matrixsize()) {
        continue;
      }
      for (Index k = 0; k < n; ++k) {
        box.AddtoBigMatrix(result[std::size_t(k)],
                           Block(i, k, box.Matrixsize()));
      }
    }
    return result;
  }
  std::vector<Eigen::MatrixXd*> targets;
  for (Eigen::MatrixXd& m : result) {
    targets.push_back(&m);
  }
  GridBoxAccumulator acc(targets, budget);
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    std::vector<Eigen::MatrixXd> blocks;
    for (Index k = 0; k < n; ++k) {
      // every third matrix gets no contribution from odd boxes
      blocks.push_back((k % 3 == 2 && i % 2 == 1)
                           ? Eigen::MatrixXd()
                           : Block(i, k, box.Matrixsize()));
    }
    std::vector<const Eigen::MatrixXd*> ptr;
    for (const Eigen::MatrixXd& m : blocks) {
      ptr.push_back(&m);
    }
    acc.Add(box, ptr);
  }
  acc.Finish();
  return result;
}
}  // namespace

BOOST_AUTO_TEST_CASE(threadlocal_and_locks_equal_serial) {
  libint2::initialize();
  QMMolecule mol("water", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/chargeresponse/water.xyz");
  BasisSet bs;
  bs.Load(std::string(XTP_TEST_DATA_FOLDER) + "/hirshfeldpartition/3-21G.xml");
  AOBasis basis;
  basis.Fill(bs, mol);
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, basis);
  const Index n_ao = basis.AOBasisSize();
  const Index n = 5;

  // serial reference, with the empty blocks skipped as in Accumulate
  std::vector<Eigen::MatrixXd> ref(std::size_t(n),
                                   Eigen::MatrixXd::Zero(n_ao, n_ao));
  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    for (Index k = 0; k < n; ++k) {
      if (k % 3 == 2 && i % 2 == 1) {
        continue;
      }
      box.AddtoBigMatrix(ref[std::size_t(k)], Block(i, k, box.Matrixsize()));
    }
  }
  const std::vector<Eigen::MatrixXd> local = Accumulate(grid, n_ao, n, 1e12);
  const std::vector<Eigen::MatrixXd> locked = Accumulate(grid, n_ao, n, 0.0);
  for (Index k = 0; k < n; ++k) {
    const double scale = ref[std::size_t(k)].cwiseAbs().maxCoeff();
    BOOST_CHECK(scale > 0.0);
    BOOST_CHECK_SMALL(
        (local[std::size_t(k)] - ref[std::size_t(k)]).cwiseAbs().maxCoeff() /
            scale,
        1e-12);
    BOOST_CHECK_SMALL(
        (locked[std::size_t(k)] - ref[std::size_t(k)]).cwiseAbs().maxCoeff() /
            scale,
        1e-12);
  }
  // the strategy follows the budget (with more than one thread)
  std::vector<Eigen::MatrixXd> m(1, Eigen::MatrixXd::Zero(n_ao, n_ao));
  GridBoxAccumulator big({&m[0]}, 1e12);
  GridBoxAccumulator none({&m[0]}, 0.0);
  BOOST_CHECK(!none.UsesThreadLocalCopies());
  BOOST_CHECK(big.UsesThreadLocalCopies() == (OPENMP::getMaxThreads() > 1));
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
