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
#ifndef VOTCA_XTP_BASISSET_H
#define VOTCA_XTP_BASISSET_H

// Standard includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace xtp {

enum class L { S = 0, P = 1, D = 2, F = 3, G = 4, H = 5, I = 6 };

std::string EnumToString(L l);

L StringToEnum(const std::string& type);
L StringToEnum(char type);

// shell type (S, P, D))

bool CheckShellType(const std::string& shelltype);

Index OffsetFuncShell(L l);

Index NumFuncShell(L l);
Index NumFuncShell_cartesian(L l);

Index OffsetFuncShell_cartesian(L l);

// Gaussian function: contraction*exp(-decay*r^2)
class GaussianPrimitive {
 public:
  GaussianPrimitive(double decay, double contraction)
      : _decay(decay), _contraction(contraction) {}
  double contraction() const { return _contraction; }

  double decay() const { return _decay; }

 private:
  double _decay;
  double _contraction;
};

class Shell {

 public:
  Shell(L l, double scale) : _l(l), _scale(scale) { ; }

  L getL() const { return _l; }

  Index getnumofFunc() const { return NumFuncShell(_l); };

  Index getOffset() const { return OffsetFuncShell(_l); }

  double getScale() const { return _scale; }

  Index getSize() const { return _gaussians.size(); }

  std::vector<GaussianPrimitive>::const_iterator begin() const {
    return _gaussians.begin();
  }
  std::vector<GaussianPrimitive>::const_iterator end() const {
    return _gaussians.end();
  }

  // adds a Gaussian
  GaussianPrimitive& addGaussian(double decay, double contraction);
  friend std::ostream& operator<<(std::ostream& out, const Shell& shell);

 private:
  L _l;
  // scaling factor
  double _scale;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<GaussianPrimitive> _gaussians;
};

/*
 * A collection of shells associated with a specific element
 */
class Element {

 public:
  Element(std::string type) : _type(type) { ; }
  using ShellIterator = std::vector<Shell>::const_iterator;
  ShellIterator begin() const { return _shells.begin(); }
  ShellIterator end() const { return _shells.end(); }

  const std::string& getType() const { return _type; }

  Shell& addShell(L l, double shellScale) {
    _shells.push_back(Shell(l, shellScale));
    return _shells.back();
  }

  Index NumOfShells() const { return _shells.size(); }

  friend std::ostream& operator<<(std::ostream& out, const Element& element);

 private:
  std::string _type;
  std::vector<Shell> _shells;
};

/*
 * A collection of elements and shells forms the basis set
 */
class BasisSet {
 public:
  void Load(const std::string& name);

  const Element& getElement(std::string element_type) const;

  std::map<std::string, Element>::iterator begin() { return _elements.begin(); }
  std::map<std::string, Element>::iterator end() { return _elements.end(); }

  std::map<std::string, Element>::const_iterator begin() const {
    return _elements.begin();
  }
  std::map<std::string, Element>::const_iterator end() const {
    return _elements.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const BasisSet& basis);

  const std::string& Name() const { return _name; }

 private:
  Element& addElement(std::string elementType);
  std::string _name;
  std::map<std::string, Element> _elements;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BASISSET_H
