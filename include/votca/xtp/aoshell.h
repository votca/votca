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

#include <libint2.hpp>

namespace votca {
namespace xtp {

class AOBasis;
class AOShell;

class AOGaussianPrimitive {
  friend AOShell;

 public:
  AOGaussianPrimitive(const GaussianPrimitive& gaussian,
                      const AOShell& aoshell);

  AOGaussianPrimitive(const AOGaussianPrimitive& gaussian,
                      const AOShell& aoshell);

  double getPowfactor() const { return _powfactor; }
  double getDecay() const { return _decay; }
  double getContraction() const { return _contraction; }
  const AOShell& getShell() const { return _aoshell; }

 private:
  static double CalcPowFactor(double decay) {
    return std::pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75);
  }
  double _decay;
  double _contraction;
  const AOShell& _aoshell;
  double _powfactor;  // used in evalspace to speed up DFT
};

/*
 * shells in a Gaussian-basis expansion
 */
class AOShell {
  friend class AOBasis;

 public:
  AOShell(const Shell& shell, const QMAtom& atom, Index startIndex);

  AOShell(const AOShell& shell);

  L getL() const { return _l; }
  Index getNumFunc() const { return NumFuncShell(_l); };
  Index getCartesianNumFunc() const { return NumFuncShell_cartesian(_l); };
  Index getStartIndex() const { return _startIndex; }
  Index getOffset() const { return OffsetFuncShell(_l); }
  Index getCartesianOffset() const { return OffsetFuncShell_cartesian(_l); }
  Index getAtomIndex() const { return _atomindex; }
  Index getSize() const { return _gaussians.size(); }

  libint2::Shell LibintShell() const;

  const Eigen::Vector3d& getPos() const { return _pos; }
  double getScale() const { return _scale; }

  void CalcMinDecay() {
    _mindecay = std::numeric_limits<double>::max();
    for (auto& gaussian : _gaussians) {
      _mindecay = std::min(_mindecay, gaussian.getDecay());
    }
  }

  double getMinDecay() const { return _mindecay; }

  void EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                   const Eigen::Vector3d& grid_pos) const;
  void EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                   Eigen::Block<Eigen::MatrixX3d>& AODervalues,
                   const Eigen::Vector3d& grid_pos) const;

  // iterator over pairs (decay constant; contraction coefficient)
  using GaussianIterator = std::vector<AOGaussianPrimitive>::const_iterator;
  GaussianIterator begin() const { return _gaussians.begin(); }
  GaussianIterator end() const { return _gaussians.end(); }

  // adds a Gaussian
  void addGaussian(const GaussianPrimitive& gaussian) {
    _gaussians.push_back(AOGaussianPrimitive(gaussian, *this));
    return;
  }

  void normalizeContraction();

  friend std::ostream& operator<<(std::ostream& out, const AOShell& shell);

 private:
  L _l;
  // scaling factor
  double _scale;
  // number of functions in shell
  double _mindecay;
  Index _startIndex;
  Eigen::Vector3d _pos;
  Index _atomindex;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<AOGaussianPrimitive> _gaussians;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOSHELL_H
