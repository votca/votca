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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/aopotential.h"
#include "votca/xtp/aotransform.h"
#include "votca/xtp/ecpaobasis.h"
#include <libecpint/ecpint.hpp>
#include <libecpint/gshell.hpp>
#include <libecpint/mathutil.hpp>
#include <libint2/solidharmonics.h>

namespace votca {
namespace xtp {

std::ostream& operator<<(std::ostream& out,
                         const libecpint::GaussianShell& shell) {
  out << " Shelltype:" << xtp::EnumToString(static_cast<L>(shell.am()))
      << " L:" << Index(shell.am()) << " Func:" << shell.nprimitive() << "\n";
  for (int i = 0; i < shell.nprimitive(); i++) {
    out << " Gaussian Decay: " << shell.exp(i);
    out << " Contractions: " << shell.coef(i) << "\n";
  }
  return out;
}

using MatrixLibInt =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void AOECP::FillPotential(const AOBasis& aobasis, ECPAOBasis ecp) {

  _aopotential =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  std::vector<libecpint::GaussianShell> basis;
  std::vector<Index> cartesian_size;
  std::vector<Index> spherical_size;
  for (const auto& shell : aobasis) {
    libecpint::GaussianShell s(
        {shell.getPos().x(), shell.getPos().y(), shell.getPos().z()},
        int(shell.getL()));
    spherical_size.push_back(shell.getNumFunc());
    cartesian_size.push_back(shell.getCartesianNumFunc());
    //The normalisation libecpint requires is identical to the libint normalisation of shells
    libint2::Shell s_libint=shell.LibintShell();
    for (Index i=0;i<Index(s_libint.nprim());i++) {
      s.addPrim(s_libint.alpha[i],s_libint.contr[0].coeff[i]);
    }
    basis.push_back(s);
  }
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  std::vector<libecpint::ECPIntegral> engines(
      OPENMP::getMaxThreads(),
      libecpint::ECPIntegral(int(aobasis.getMaxL()), int(ecp.getMaxL()), 0));
  #pragma omp parallel for schedule(guided)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libecpint::ECPIntegral& engine = engines[thread_id];
    Index bf1 = shell2bf[s1];
    Index n1 = spherical_size[s1];
    Index c1 = cartesian_size[s1];
    for (Index s2 = 0; s2 <= s1; ++s2) {
      Index bf2 = shell2bf[s2];
      Index n2 = spherical_size[s2];
      Index c2 = cartesian_size[s2];

      MatrixLibInt cartesian_result = MatrixLibInt::Zero(c1, c2);
      for (auto& ecppotential : ecp) {
        libecpint::TwoIndex<double> results;
        engine.compute_shell_pair(ecppotential, basis[s1], basis[s2], results);
        cartesian_result +=
            Eigen::Map<MatrixLibInt>(results.data.data(), c1, c2);
      }

      MatrixLibInt spherical_result = MatrixLibInt::Zero(n1, n2);
      libint2::solidharmonics::tform<double>(basis[s1].l, basis[s2].l,
                                             cartesian_result.data(),
                                             spherical_result.data());
      _aopotential.block(bf1, bf2, n1, n2) = spherical_result;
      if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                       // {s2,s1} block, note the transpose!
        _aopotential.block(bf2, bf1, n2, n1) = spherical_result.transpose();
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca
