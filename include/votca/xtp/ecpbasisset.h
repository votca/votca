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
#ifndef VOTCA_XTP_PSEUDOPOTENTIAL_H
#define VOTCA_XTP_PSEUDOPOTENTIAL_H

#include <votca/xtp/basisset.h>

namespace votca {
namespace xtp {

class ECPGaussianPrimitive {
 public:
  ECPGaussianPrimitive(Index power, double decay, double contraction)
      : _power(power), _decay(decay), _contraction(contraction) {
    ;
  }

  Index _power;
  double _decay;
  double _contraction;
};

class ECPShell {

 public:
  ECPShell(std::string type) : _type(type) { ; }
  const std::string& getType() const { return _type; }

  Index getL() const { return FindLmax(_type); }

  Index getnumofFunc() const { return NumFuncShell(_type); };

  Index getOffset() const { return OffsetFuncShell(_type); }

  Index getSize() const { return _gaussians.size(); }

  std::vector<ECPGaussianPrimitive>::const_iterator begin() const {
    return _gaussians.begin();
  }
  std::vector<ECPGaussianPrimitive>::const_iterator end() const {
    return _gaussians.end();
  }

  // adds a Gaussian of a pseudopotential
  ECPGaussianPrimitive& addGaussian(Index power, double decay,
                                    double contraction);

  friend std::ostream& operator<<(std::ostream& out, const ECPShell& shell);

 private:
  std::string _type;
  // vector of pairs of decay constants and contraction coefficients
  std::vector<ECPGaussianPrimitive> _gaussians;
};

/*
 * A collection of shells associated with a specific element
 */
class ECPElement {
 public:
  ECPElement(std::string type, Index lmax, Index ncore)
      : _type(type), _lmax(lmax), _ncore(ncore) {
    ;
  }
  using ECPShellIterator = std::vector<ECPShell>::const_iterator;
  ECPShellIterator begin() const { return _shells.begin(); }
  ECPShellIterator end() const { return _shells.end(); }

  const std::string& getType() const { return _type; }

  Index getLmax() const { return _lmax; }

  Index getNcore() const { return _ncore; }

  ECPShell& addShell(const std::string& shellType) {
    _shells.push_back(ECPShell(shellType));
    return _shells.back();
  }

  Index NumOfShells() const { return _shells.size(); }

  friend std::ostream& operator<<(std::ostream& out, const ECPElement& element);

 private:
  std::string _type;
  //  applies to the highest angular momentum lmax
  Index _lmax;
  // replaces ncore electrons
  Index _ncore;

  std::vector<ECPShell> _shells;
};

/*
 * A collection of elements and shells forms the basis set
 */
class ECPBasisSet {
 public:
  void Load(const std::string& name);

  ECPElement& addElement(std::string elementType, Index lmax, Index ncore);

  const ECPElement& getElement(std::string element_type) const;

  std::map<std::string, std::shared_ptr<ECPElement> >::iterator begin() {
    return _elements.begin();
  }
  std::map<std::string, std::shared_ptr<ECPElement> >::iterator end() {
    return _elements.end();
  }

  std::map<std::string, std::shared_ptr<ECPElement> >::const_iterator begin()
      const {
    return _elements.begin();
  }
  std::map<std::string, std::shared_ptr<ECPElement> >::const_iterator end()
      const {
    return _elements.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const ECPBasisSet& basis);

 private:
  std::string _name;
  std::map<std::string, std::shared_ptr<ECPElement> > _elements;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PSEUDOPOTENTIAL_H
