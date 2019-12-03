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

#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>
#include <votca/xtp/eigen.h>

#include "qmatom.h"
#include <votca/xtp/basisset.h>

namespace votca {
namespace xtp {

class AOBasis;
class AOShell;

class AOGaussianPrimitive {
  friend class AOShell;

 public:
  double getPowfactor() const { return _powfactor; }
  double getDecay() const { return _decay; }
  const Eigen::VectorXd& getContraction() const { return _contraction; }
  const AOShell& getShell() const { return _aoshell; }

 private:
  double _decay;
  Eigen::VectorXd _contraction;
  const AOShell& _aoshell;
  double _powfactor;  // used in evalspace to speed up DFT
  // private constructor, only a shell can create a primitive
  AOGaussianPrimitive(const GaussianPrimitive& gaussian, const AOShell& aoshell)
      : _decay(gaussian.decay()), _aoshell(aoshell) {
    _contraction = Eigen::VectorXd::Map(gaussian.Contractions().data(),
                                        gaussian.Contractions().size());
    _powfactor =
        std::pow(2.0 * _decay / boost::math::constants::pi<double>(), 0.75);
  }

  AOGaussianPrimitive(const AOGaussianPrimitive& gaussian,
                      const AOShell& aoshell)
      : _decay(gaussian._decay),
        _contraction(gaussian._contraction),
        _aoshell(aoshell),
        _powfactor(gaussian._powfactor) {
    ;
  }
};

/*
 * shells in a Gaussian-basis expansion
 */
class AOShell {
  friend class AOBasis;

 public:
  AOShell(const AOShell& shell) {

    _type = shell._type;
    _Lmax = shell._Lmax;
    _Lmin = shell._Lmin;
    _scale = shell._scale;
    _numFunc = shell._numFunc;
    _numcartFunc = shell._numcartFunc;
    _mindecay = shell._mindecay;
    _startIndex = shell._startIndex;
    _offset = shell._offset;
    _cartOffset = shell._cartOffset;
    _pos = shell._pos;
    _atomindex = shell._atomindex;
    _gaussians.reserve(shell._gaussians.size());
    for (const auto& gaus : shell._gaussians) {
      _gaussians.push_back(AOGaussianPrimitive(gaus, *this));
    }
  }

  const std::string& getType() const { return _type; }
  Index getNumFunc() const { return _numFunc; }
  Index getCartesianNumFunc() const { return _numcartFunc; }
  Index getStartIndex() const { return _startIndex; }
  Index getOffset() const { return _offset; }
  Index getCartesianOffset() const { return _cartOffset; }
  Index getAtomIndex() const { return _atomindex; }

  Index getLmax() const { return _Lmax; }
  Index getLmin() const { return _Lmin; }

  bool isCombined() const { return _Lmax != _Lmin; }

  const Eigen::Vector3d& getPos() const { return _pos; }
  double getScale() const { return _scale; }

  void CalcMinDecay() {
    _mindecay = std::numeric_limits<double>::max();
    for (auto& gaussian : _gaussians) {
      if (gaussian.getDecay() < _mindecay) {
        _mindecay = gaussian.getDecay();
      }
    }
    return;
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
  // only class aobasis can construct shells
  AOShell(const Shell& shell, const QMAtom& atom, Index startIndex)
      : _type(shell.getType()),
        _Lmax(shell.getLmax()),
        _Lmin(shell.getLmin()),
        _scale(shell.getScale()),
        _numFunc(shell.getnumofFunc()),
        _numcartFunc(xtp::NumFuncShell_cartesian(shell.getType())),
        _startIndex(startIndex),
        _offset(shell.getOffset()),
        _cartOffset(xtp::OffsetFuncShell_cartesian(shell.getType())),
        _pos(atom.getPos()),
        _atomindex(atom.getId()) {
    ;
  }

  std::string _type;
  Index _Lmax;
  Index _Lmin;
  // scaling factor
  double _scale;
  // number of functions in shell
  Index _numFunc;
  Index _numcartFunc;
  double _mindecay;
  Index _startIndex;
  Index _offset;
  Index _cartOffset;
  Eigen::Vector3d _pos;
  Index _atomindex;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<AOGaussianPrimitive> _gaussians;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOSHELL_H
