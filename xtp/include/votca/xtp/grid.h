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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_GRID_H
#define VOTCA_XTP_GRID_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "logger.h"
#include "qmatom.h"
#include "qmmolecule.h"
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
  const std::vector<Eigen::Vector3d> &getGridPositions() const {
    return gridpoints_;
  }

  Eigen::VectorXd &getGridValues() { return gridvalues_; }
  const Eigen::VectorXd &getGridValues() const { return gridvalues_; }
  Index size() { return Index(gridpoints_.size()); }

  void printGridtoxyzfile(std::string filename);

  void setupCHELPGGrid(const QMMolecule &Atomlist) {
    padding_ = 3 * tools::conv::ang2bohr;  // Additional distance from molecule
                                           // to set up grid according to CHELPG
                                           // paper [Journal of Computational
                                           // Chemistry 11, 361, 1990]
    gridspacing_ =
        0.3 * tools::conv::ang2bohr;  // Grid spacing according to same paper
    cutoff_ = 2.8 * tools::conv::ang2bohr;
    setupgrid(Atomlist);
  }

 private:
  void setupgrid(const QMMolecule &Atomlist);
  std::vector<Eigen::Vector3d> gridpoints_;
  Eigen::VectorXd gridvalues_;

  double cutoff_;
  double gridspacing_;
  double padding_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GRID_H
