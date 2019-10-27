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
#ifndef VOTCA_XTP_ECPAOSHELL_H
#define VOTCA_XTP_ECPAOSHELL_H

#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>
#include <votca/xtp/eigen.h>

#include "qmatom.h"
#include <votca/xtp/ecpbasisset.h>

namespace votca {
namespace xtp {

class ECPAOGaussianPrimitive {

 public:
  ECPAOGaussianPrimitive(const ECPGaussianPrimitive& gaussian)
      : _power(gaussian._power),
        _decay(gaussian._decay),
        _contraction(gaussian._contraction) {
    ;
  }
  int getPower() const { return _power; }
  double getDecay() const { return _decay; }
  double getContraction() const { return _contraction; }

 private:
  int _power = 0;
  double _decay = 0.0;
  double _contraction = 0.0;
};

/*
 * shells in a Gaussian-basis expansion
 */
class ECPAOShell {
 public:
  ECPAOShell(const ECPShell& shell, const QMAtom& atom, int startIndex,
             int Lmax)
      : _type(shell.getType()),
        _L(shell.getL()),
        _numFunc(shell.getnumofFunc()),
        _startIndex(startIndex),
        _offset(shell.getOffset()),
        _pos(atom.getPos()),
        _atomindex(atom.getId()),
        _Lmax_element(Lmax) {
    ;
  }

  const std::string& getType() const { return _type; }
  int getNumFunc() const { return _numFunc; }
  int getStartIndex() const { return _startIndex; }
  int getOffset() const { return _offset; }
  long getAtomIndex() const { return _atomindex; }

  int getL() const { return _L; }
  int getLmaxElement() const { return _Lmax_element; }
  // Local part is with L=Lmax
  bool isNonLocal() const { return (_L < _Lmax_element); }
  const Eigen::Vector3d& getPos() const { return _pos; }

  long getSize() const { return _gaussians.size(); }

  // iterator over pairs (decay constant; contraction coefficient)
  using ECPGaussianIterator =
      std::vector<ECPAOGaussianPrimitive>::const_iterator;
  ECPGaussianIterator begin() const { return _gaussians.begin(); }
  ECPGaussianIterator end() const { return _gaussians.end(); }

  // adds a Gaussian
  void addGaussian(const ECPGaussianPrimitive& gaussian) {
    _gaussians.push_back(ECPAOGaussianPrimitive(gaussian));
    return;
  }

  friend std::ostream& operator<<(std::ostream& out, const ECPAOShell& shell);

 private:
  std::string _type;
  int _L;
  // number of functions in shell
  int _numFunc;
  int _startIndex;
  int _offset;
  Eigen::Vector3d _pos;
  long _atomindex;
  int _Lmax_element;  // Lmax of the Element not the shell

  // vector of pairs of decay constants and contraction coefficients
  std::vector<ECPAOGaussianPrimitive> _gaussians;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ECPAOSHELL_H
