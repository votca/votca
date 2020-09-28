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
#include <vector>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {

std::unordered_map<Index, std::vector<Index>> AOMatrix::compute_shellpairs(
    const AOBasis& aobasis, const double threshold) {

  Index nthreads = OPENMP::getMaxThreads();

  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  Index nsh1 = shells.size();

  // construct the 2-electron repulsion integrals engine
  std::vector<libint2::Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(libint2::Operator::overlap, aobasis.getMaxNprim(),
                       aobasis.getMaxL(), 0);
  for (Index i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::unordered_map<Index, std::vector<Index>> splist;

  libint2::Engine& engine = engines[thread_id];
  const libint2::Engine::target_ptr_vec& buf = engine.results();

  // loop over permutationally-unique set of shells
  for (Index s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {

#pragma omp critical
    {
      if (splist.find(s1) == splist.end())
        splist.insert(std::make_pair(s1, std::vector<Index>()));
    }

    Index n1 = shells[s1].size();  // number of basis functions in this shell

    for (Index s2 = 0; s2 <= s1; ++s2, ++s12) {
      if (s12 % nthreads != thread_id) continue;

      bool on_same_center = (shells[s1].O == shells[s2].O);
      bool significant = on_same_center;
      if (not on_same_center) {
        Index n2 = shells[s2].size();
        engines[thread_id].compute(shells[s1], shells[s2]);
        Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
        double norm = buf_mat.norm();
        significant = (norm >= threshold);
      }

      if (significant) {
#pragma omp critical
        { splist[s1].emplace_back(s2); }
      }
    }
  }
};  // end of compute

parallel_do(compute);

// resort shell list in increasing order, i.e. splist[s][s1] < splist[s][s2]
// if s1 < s2 N.B. only parallelized over 1 shell index
auto sort = [&](Index thread_id) {
  for (Index s1 = 0l; s1 != nsh1; ++s1) {
    if (s1 % nthreads == thread_id) {
      std::vector<Index>& list = splist[s1];
      std::sort(list.begin(), list.end());
    }
  }
};  // end of sort

parallel_do(sort);

return splist;
}  // namespace xtp

}  // namespace votca
}  // namespace votca
