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

#ifndef __XTP_AOSHELL__H
#define __XTP_AOSHELL__H

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
  int getPower() const { return _power; }
  double getDecay() const { return _decay; }
  const std::vector<double>& getContraction() const { return _contraction; }
  const AOShell* getShell() const { return _aoshell; }

 private:
  int _power;  // used in pseudopotenials only
  double _decay;
  std::vector<double> _contraction;
  AOShell* _aoshell;
  double _powfactor;  // used in evalspace to speed up DFT
  // private constructor, only a shell can create a primitive
  AOGaussianPrimitive(const GaussianPrimitive& gaussian, AOShell* aoshell)
      : _power(gaussian._power),
        _decay(gaussian._decay),
        _contraction(gaussian._contraction),
        _aoshell(aoshell) {
    _powfactor =
        std::pow(2.0 * _decay / boost::math::constants::pi<double>(), 0.75);
  }

  AOGaussianPrimitive(const AOGaussianPrimitive& gaussian, AOShell* aoshell)
      : _power(gaussian._power),
        _decay(gaussian._decay),
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
    _scale = shell._scale;
    _numFunc = shell._numFunc;
    _mindecay = shell._mindecay;
    _startIndex = shell._startIndex;
    _offset = shell._offset;
    _pos = shell._pos;
    _atomindex = shell._atomindex;
    _nonlocal = shell._nonlocal;
    _gaussians.reserve(shell._gaussians.size());
    for (const auto& gaus : shell._gaussians) {
      _gaussians.push_back(AOGaussianPrimitive(gaus, this));
    }
  }

  const std::string& getType() const { return _type; }
  int getNumFunc() const { return _numFunc; }
  int getStartIndex() const { return _startIndex; }
  int getOffset() const { return _offset; }
  int getAtomIndex() const { return _atomindex; }

  int getLmax() const { return _Lmax; }

  bool isCombined() const { return _type.length() > 1; }

  bool isNonLocal() const { return _nonlocal; }

  const tools::vec& getPos() const { return _pos; }
  double getScale() const { return _scale; }

  int getSize() const { return _gaussians.size(); }

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
                   const tools::vec& grid_pos) const;
  void EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                   Eigen::Block<Eigen::MatrixX3d>& AODervalues,
                   const tools::vec& grid_pos) const;

  // iterator over pairs (decay constant; contraction coefficient)
  typedef std::vector<AOGaussianPrimitive>::const_iterator GaussianIterator;
  GaussianIterator begin() const { return _gaussians.begin(); }
  GaussianIterator end() const { return _gaussians.end(); }

  // adds a Gaussian
  void addGaussian(const GaussianPrimitive& gaussian) {
    AOGaussianPrimitive aogaussian = AOGaussianPrimitive(gaussian, this);
    _gaussians.push_back(aogaussian);
    return;
  }

  void normalizeContraction();

  friend std::ostream& operator<<(std::ostream& out, const AOShell& shell);

 private:
  // only class aobasis can construct shells
  AOShell(const Shell& shell, const QMAtom& atom, int startIndex)
      : _type(shell.getType()),
        _Lmax(shell.getLmax()),
        _scale(shell.getScale()),
        _numFunc(shell.getnumofFunc()),
        _startIndex(startIndex),
        _offset(shell.getOffset()),
        _pos(atom.getPos()),
        _atomindex(atom.getAtomID()) {
    ;
  }
  // for ECPs
  AOShell(const Shell& shell, const QMAtom& atom, int startIndex, bool nonlocal)
      : _type(shell.getType()),
        _Lmax(shell.getLmax()),
        _scale(shell.getScale()),
        _numFunc(shell.getnumofFunc()),
        _startIndex(startIndex),
        _offset(shell.getOffset()),
        _pos(atom.getPos()),
        _atomindex(atom.getAtomID()),
        _nonlocal(nonlocal) {
    ;
  }

  // only class aobasis can destruct shells
  ~AOShell(){};

  // shell type (S, P, D))
  std::string _type;
  int _Lmax;
  // scaling factor
  double _scale;
  // number of functions in shell
  int _numFunc;
  double _mindecay;
  int _startIndex;
  int _offset;
  tools::vec _pos;
  int _atomindex;
  // used for ecp calculations
  bool _nonlocal;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<AOGaussianPrimitive> _gaussians;
};

}  // namespace xtp
}  // namespace votca

#endif /* AOSHELL_H */
