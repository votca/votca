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

#ifndef VOTCA_XTP_BASISSET_H
#define VOTCA_XTP_BASISSET_H

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace votca {
namespace xtp {
// shell type (S, P, D))

int FindLmax(const std::string& type);

int FindLmin(const std::string& type);

int OffsetFuncShell(const std::string& shell_type);

int NumFuncShell(const std::string& shell_type);
int NumFuncShell_cartesian(const std::string& shell_type);

int OffsetFuncShell_cartesian(const std::string& shell_type);

std::vector<int> NumFuncSubShell(const std::string& shell_type);

class Shell;
class Element;
class BasisSet;

// Gaussian function: contraction*exp(-decay*r^2)
class GaussianPrimitive {
  friend class Shell;

 public:
  int _power;  // used in pseudopotenials only
  double _decay;
  std::vector<double> _contraction;

 private:
  // private constructor, only a shell can create a primitive
  GaussianPrimitive(double decay, std::vector<double> contraction)
      : _decay(decay), _contraction(contraction) {
    ;
  }

  GaussianPrimitive(int power, double decay, std::vector<double> contraction)
      : _power(power), _decay(decay), _contraction(contraction) {
    ;
  }
};

class Shell {
  friend class Element;

 public:
  const std::string& getType() const { return _type; }

  bool isCombined() const { return (_type.length() > 1); }

  int getLmax() const { return FindLmax(_type); }

  int getnumofFunc() const { return NumFuncShell(_type); };

  int getOffset() const { return OffsetFuncShell(_type); }

  double getScale() const { return _scale; }

  int getSize() const { return _gaussians.size(); }

  std::vector<GaussianPrimitive>::const_iterator begin() const {
    return _gaussians.begin();
  }
  std::vector<GaussianPrimitive>::const_iterator end() const {
    return _gaussians.end();
  }

  // adds a Gaussian
  GaussianPrimitive& addGaussian(double decay, std::vector<double> contraction);

  // adds a Gaussian of a pseudopotential
  GaussianPrimitive& addGaussian(int power, double decay,
                                 std::vector<double> contraction);

  friend std::ostream& operator<<(std::ostream& out, const Shell& shell);

 private:
  // only class Element can construct shells
  Shell(std::string type, double scale) : _type(type), _scale(scale) { ; }

  std::string _type;
  // scaling factor
  double _scale;

  // vector of pairs of decay constants and contraction coefficients
  std::vector<GaussianPrimitive> _gaussians;
};

/*
 * A collection of shells associated with a specific element
 */
class Element {
  friend class BasisSet;

 public:
  typedef std::vector<Shell>::const_iterator ShellIterator;
  ShellIterator begin() const { return _shells.begin(); }
  ShellIterator end() const { return _shells.end(); }

  const std::string& getType() const { return _type; }

  int getLmax() const { return _lmax; }

  int getNcore() const { return _ncore; }

  Shell& addShell(const std::string& shellType, double shellScale) {
    _shells.push_back(Shell(shellType, shellScale));
    return _shells.back();
  }

  int NumOfShells() const { return _shells.size(); }

  friend std::ostream& operator<<(std::ostream& out, const Element& element);

 private:
  // only class BasisSet can create Elements
  Element(std::string type) : _type(type) { ; }

  // used for the pseudopotential
  Element(std::string type, int lmax, int ncore)
      : _type(type), _lmax(lmax), _ncore(ncore) {
    ;
  }

  // only class BasisSet can destruct Elements

  std::string _type;
  // lmax is used in the pseudopotentials only (applies to the highest angular
  // momentum lmax)
  int _lmax;
  // ncore is used in the pseudopotentials only (replaces ncore electrons))
  int _ncore;

  std::vector<Shell> _shells;
};

/*
 * A collection of elements and shells forms the basis set
 */
class BasisSet {
 public:
  void LoadBasisSet(const std::string& name);

  void LoadPseudopotentialSet(const std::string& name);

  Element& addElement(std::string elementType);

  // used for pseudopotentials only
  Element& addElement(std::string elementType, int lmax, int ncore);

  const Element& getElement(std::string element_type) const;

  std::map<std::string, std::shared_ptr<Element> >::iterator begin() {
    return _elements.begin();
  }
  std::map<std::string, std::shared_ptr<Element> >::iterator end() {
    return _elements.end();
  }

  std::map<std::string, std::shared_ptr<Element> >::const_iterator begin()
      const {
    return _elements.begin();
  }
  std::map<std::string, std::shared_ptr<Element> >::const_iterator end() const {
    return _elements.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const BasisSet& basis);

 private:
  std::string _name;
  std::map<std::string, std::shared_ptr<Element> > _elements;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BASISSET_H
