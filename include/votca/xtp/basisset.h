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
      : decay_(decay), contraction_(contraction) {}
  double contraction() const { return contraction_; }

  double decay() const { return decay_; }

 private:
  double decay_;
  double contraction_;
};

class Shell {

 public:
  Shell(L l, double scale) : l_(l), scale_(scale) { ; }

  L getL() const { return l_; }

  Index getnumofFunc() const { return NumFuncShell(l_); };

  Index getOffset() const { return OffsetFuncShell(l_); }

  double getScale() const { return scale_; }

  Index getSize() const { return gaussians_.size(); }

  std::vector<GaussianPrimitive>::const_iterator begin() const {
    return gaussians_.begin();
  }
  std::vector<GaussianPrimitive>::const_iterator end() const {
    return gaussians_.end();
  }

  // adds a Gaussian
  GaussianPrimitive& addGaussian(double decay, double contraction);
  friend std::ostream& operator<<(std::ostream& out, const Shell& shell);

 private:
  L l_;
  // scaling factor
  double scale_;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<GaussianPrimitive> gaussians_;
};

/*
 * A collection of shells associated with a specific element
 */
class Element {

 public:
  Element(std::string type) : type_(type) { ; }
  using ShellIterator = std::vector<Shell>::const_iterator;
  ShellIterator begin() const { return shells_.begin(); }
  ShellIterator end() const { return shells_.end(); }

  const std::string& getType() const { return type_; }

  Shell& addShell(L l, double shellScale) {
    shells_.push_back(Shell(l, shellScale));
    return shells_.back();
  }

  Index NumOfShells() const { return shells_.size(); }

  friend std::ostream& operator<<(std::ostream& out, const Element& element);

 private:
  std::string type_;
  std::vector<Shell> shells_;
};

/*
 * A collection of elements and shells forms the basis set
 */
class BasisSet {
 public:
  void Load(const std::string& name);

  const Element& getElement(std::string element_type) const;

  std::map<std::string, Element>::iterator begin() { return elements_.begin(); }
  std::map<std::string, Element>::iterator end() { return elements_.end(); }

  std::map<std::string, Element>::const_iterator begin() const {
    return elements_.begin();
  }
  std::map<std::string, Element>::const_iterator end() const {
    return elements_.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const BasisSet& basis);

  const std::string& Name() const { return name_; }

 private:
  Element& addElement(std::string elementType);
  std::string name_;
  std::map<std::string, Element> elements_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BASISSET_H
