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
#ifndef VOTCA_XTP_ECPBASISSET_H
#define VOTCA_XTP_ECPBASISSET_H

// Local VOTCA includes
#include "basisset.h"

namespace votca {
namespace xtp {

class ECPGaussianPrimitive {
 public:
  ECPGaussianPrimitive(Index power, double decay, double contraction)
      : power_(power), decay_(decay), contraction_(contraction) {
    ;
  }

  Index power_;
  double decay_;
  double contraction_;
};

class ECPShell {

 public:
  ECPShell(L l) : l_(l) { ; }
  L getL() const { return l_; }

  Index getnumofFunc() const { return NumFuncShell(l_); }

  Index getOffset() const { return OffsetFuncShell(l_); }

  Index getSize() const { return gaussians_.size(); }

  std::vector<ECPGaussianPrimitive>::const_iterator begin() const {
    return gaussians_.begin();
  }
  std::vector<ECPGaussianPrimitive>::const_iterator end() const {
    return gaussians_.end();
  }

  // adds a Gaussian of a pseudopotential
  ECPGaussianPrimitive& addGaussian(Index power, double decay,
                                    double contraction);

  friend std::ostream& operator<<(std::ostream& out, const ECPShell& shell);

 private:
  L l_;
  // vector of pairs of decay constants and contraction coefficients
  std::vector<ECPGaussianPrimitive> gaussians_;
};

/*
 * A collection of shells associated with a specific element
 */
class ECPElement {
 public:
  ECPElement(std::string type, L lmax, Index ncore)
      : type_(type), lmax_(lmax), ncore_(ncore) {
    ;
  }
  using ECPShellIterator = std::vector<ECPShell>::const_iterator;
  ECPShellIterator begin() const { return shells_.begin(); }
  ECPShellIterator end() const { return shells_.end(); }

  const std::string& getType() const { return type_; }

  L getLmax() const { return lmax_; }

  Index getNcore() const { return ncore_; }

  ECPShell& addShell(L l) {
    shells_.push_back(ECPShell(l));
    return shells_.back();
  }

  Index NumOfShells() const { return shells_.size(); }

  friend std::ostream& operator<<(std::ostream& out, const ECPElement& element);

 private:
  std::string type_;
  //  applies to the highest angular momentum lmax
  L lmax_;
  // replaces ncore electrons
  Index ncore_;

  std::vector<ECPShell> shells_;
};

/*
 * A collection of elements and shells forms the basis set
 */
class ECPBasisSet {
 public:
  void Load(const std::string& name);

  ECPElement& addElement(std::string elementType, L lmax, Index ncore);

  const ECPElement& getElement(std::string element_type) const;

  std::map<std::string, std::shared_ptr<ECPElement> >::iterator begin() {
    return elements_.begin();
  }
  std::map<std::string, std::shared_ptr<ECPElement> >::iterator end() {
    return elements_.end();
  }

  const std::string& Name() const { return name_; }

  std::map<std::string, std::shared_ptr<ECPElement> >::const_iterator begin()
      const {
    return elements_.begin();
  }
  std::map<std::string, std::shared_ptr<ECPElement> >::const_iterator end()
      const {
    return elements_.end();
  }

  friend std::ostream& operator<<(std::ostream& out, const ECPBasisSet& basis);

 private:
  std::string name_;
  std::map<std::string, std::shared_ptr<ECPElement> > elements_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ECPBASISSET_H
