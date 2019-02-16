/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __XTP_GRID__H
#define __XTP_GRID__H

#include <string>
#include <vector>
#include <votca/ctp/apolarsite.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/polarseg.h>
#include <votca/tools/elements.h>
#include <votca/xtp/qmatom.h>
/**
 * \brief Takes a list of atoms, and creates CHELPG grid.
 *
 *
 *
 */

namespace votca {
namespace xtp {

class Grid {
 public:
  const std::vector<tools::vec>& getGridPositions() const {
    return _gridpoints;
  }

  Eigen::VectorXd& getGridValues() { return _gridvalues; }
  const Eigen::VectorXd& getGridValues() const { return _gridvalues; }

  unsigned getsize() { return _gridpoints.size(); }

  void printGridtoxyzfile(std::string filename);

  void setupCHELPGGrid(std::vector<QMAtom*>& Atomlist) {
    _padding = 3 * tools::conv::ang2bohr;  // Additional distance from molecule
                                           // to set up grid according to CHELPG
                                           // paper [Journal of Computational
                                           // Chemistry 11, 361, 1990]
    _gridspacing =
        0.3 * tools::conv::ang2bohr;  // Grid spacing according to same paper
    _cutoff = 2.8 * tools::conv::ang2bohr;
    setupgrid(Atomlist);
  }

 private:
  void setupgrid(std::vector<QMAtom*>& Atomlist);
  std::vector<tools::vec> _gridpoints;
  Eigen::VectorXd _gridvalues;

  double _cutoff;
  double _gridspacing;
  double _padding;
};

}  // namespace xtp
}  // namespace votca

#endif /* GRID_H */