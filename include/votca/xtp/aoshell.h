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

#pragma once
#ifndef VOTCA_XTP_AOSHELL_H
#define VOTCA_XTP_AOSHELL_H

// Third party includes
#include <boost/math/constants/constants.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "basisset.h"
#include "eigen.h"
#include "qmatom.h"
// include libint last otherwise it overrides eigen
#include <libint2/shell.h>

namespace votca {
namespace xtp {

class AOBasis;
class AOShell;
class SetupCptTable;

class AOGaussianPrimitive {

 public:
  AOGaussianPrimitive(const GaussianPrimitive& gaussian);
  friend class AOShell;

  struct data {
    Index atomid;
    Index l;
    Index startindex;
    double decay;
    double contraction;
    double x;
    double y;
    double z;
    double scale;
  };

  AOGaussianPrimitive(const AOGaussianPrimitive::data& d) {
    decay_ = d.decay;
    contraction_ = d.contraction;
    powfactor_ = CalcPowFactor(decay_);
  }

  static void SetupCptTable(CptTable& table);

  void WriteData(data& d, const AOShell& s) const;

  double getPowfactor() const { return powfactor_; }
  double getDecay() const { return decay_; }
  double getContraction() const { return contraction_; }

 private:
  static double CalcPowFactor(double decay) {
    return std::pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75);
  }
  double decay_;
  double contraction_;
  double powfactor_;  // used in evalspace to speed up DFT
};

/*
 * shells in a Gaussian-basis expansion
 */
class AOShell {
  friend AOBasis;

 public:
  AOShell(const Shell& shell, const QMAtom& atom, Index startIndex);

  AOShell(const AOGaussianPrimitive::data& d) {
    l_ = static_cast<L>(d.l);
    startIndex_ = d.startindex;
    atomindex_ = d.atomid;
    pos_ = Eigen::Vector3d(d.x, d.y, d.z);
    gaussians_.push_back(AOGaussianPrimitive(d));
  }

  L getL() const { return l_; }
  Index getNumFunc() const { return NumFuncShell(l_); };
  Index getCartesianNumFunc() const { return NumFuncShell_cartesian(l_); };
  Index getStartIndex() const { return startIndex_; }
  Index getOffset() const { return OffsetFuncShell(l_); }
  Index getCartesianOffset() const { return OffsetFuncShell_cartesian(l_); }
  Index getAtomIndex() const { return atomindex_; }
  Index getSize() const { return gaussians_.size(); }

  libint2::Shell LibintShell() const;

  const Eigen::Vector3d& getPos() const { return pos_; }

  void CalcMinDecay() {
    mindecay_ = std::numeric_limits<double>::max();
    for (auto& gaussian : gaussians_) {
      mindecay_ = std::min(mindecay_, gaussian.getDecay());
    }
  }

  double getMinDecay() const { return mindecay_; }

  struct AOValues
  {

    AOValues(Index size){
      values=Eigen::VectorXd::Zero(size);
      derivatives=Eigen::MatrixX3d::Zero(size,3);
    }
    Eigen::VectorXd values;
    Eigen::MatrixX3d derivatives;

  };

  AOValues EvalAOspace(const Eigen::Vector3d& grid_pos) const;

  // iterator over pairs (decay constant; contraction coefficient)
  using GaussianIterator = std::vector<AOGaussianPrimitive>::const_iterator;
  GaussianIterator begin() const { return gaussians_.begin(); }
  GaussianIterator end() const { return gaussians_.end(); }

  // adds a Gaussian
  void addGaussian(const GaussianPrimitive& gaussian) {
    gaussians_.push_back(AOGaussianPrimitive(gaussian));
    return;
  }

  void normalizeContraction();

  friend std::ostream& operator<<(std::ostream& out, const AOShell& shell);

 private:
  L l_;
  // scaling factor
  // number of functions in shell
  double mindecay_;
  Index startIndex_;
  Eigen::Vector3d pos_;
  Index atomindex_;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<AOGaussianPrimitive> gaussians_;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOSHELL_H
