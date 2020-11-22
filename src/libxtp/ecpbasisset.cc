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

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/tools/globals.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/ecpbasisset.h"

namespace votca {
namespace xtp {

void ECPBasisSet::Load(const std::string& name) {
  tools::Property basis_property;

  // if name contains .xml, assume a ecp .xml file is located in the working
  // directory
  std::size_t found_xml = name.find(".xml");
  std::string xmlFile;
  if (found_xml != std::string::npos) {
    xmlFile = name;
  } else {
    xmlFile = tools::GetVotcaShare() + "/xtp/ecps/" + name + ".xml";
  }
  basis_property.LoadFromXML(xmlFile);
  _name =
      basis_property.get("pseudopotential").getAttribute<std::string>("name");
  std::vector<tools::Property*> elementProps =
      basis_property.Select("pseudopotential.element");

  for (tools::Property* elementProp : elementProps) {
    std::string elementName = elementProp->getAttribute<std::string>("name");
    Index lmax = elementProp->getAttribute<Index>("lmax");
    if (lmax > Index(L::I)) {
      throw std::runtime_error("In ecps lmax larger " +
                               std::to_string(Index(L::I)) + " is not allowed");
    }
    Index ncore = elementProp->getAttribute<Index>("ncore");

    ECPElement& element = addElement(elementName, static_cast<L>(lmax), ncore);

    std::vector<tools::Property*> shellProps = elementProp->Select("shell");
    for (tools::Property* shellProp : shellProps) {
      std::string shellType = shellProp->getAttribute<std::string>("type");
      if (shellType.size() > 1) {
        throw std::runtime_error(
            "In ecps no combined shells e.g. SP are allowed. Here:" +
            shellType);
      }
      ECPShell& shell = element.addShell(StringToEnum(shellType));
      std::vector<tools::Property*> constProps = shellProp->Select("constant");
      for (tools::Property* constProp : constProps) {
        Index power = constProp->getAttribute<Index>("power");
        double decay = constProp->getAttribute<double>("decay");
        double contraction = constProp->getAttribute<double>("contraction");
        shell.addGaussian(power, decay, contraction);
      }
    }
  }
  return;
}

// adding an Element to a Pseudopotential Library
ECPElement& ECPBasisSet::addElement(std::string elementType, L lmax,
                                    Index ncore) {
  std::shared_ptr<ECPElement> element(new ECPElement(elementType, lmax, ncore));
  _elements[elementType] = element;
  return *element;
}

const ECPElement& ECPBasisSet::getElement(std::string element_type) const {
  std::map<std::string, std::shared_ptr<ECPElement> >::const_iterator itm =
      _elements.find(element_type);
  if (itm == _elements.end()) {
    throw std::runtime_error("Basis set " + _name +
                             " does not have element of type " + element_type);
  }
  const ECPElement& element = *((*itm).second);
  return element;
}

std::ostream& operator<<(std::ostream& out, const ECPShell& shell) {

  out << "Type:" << xtp::EnumToString(shell.getL())
      << " Func: " << shell.getnumofFunc() << "\n";
  for (const auto& gaussian : shell._gaussians) {
    out << " Gaussian Decay: " << gaussian._decay;
     out << " Power: " << gaussian._power;
    out << " Contraction:" << gaussian._contraction << "\n";
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const ECPElement& element) {
  out << "Element:" << element.getType() << " Lmax:" << element.getLmax()
      << " Ncore:" << element.getNcore() << "\n";
  for (const auto& shell : element) {
    out << shell;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const ECPBasisSet& basis) {
  out << "BasisSet:" << basis._name << "\n";
  for (const auto& element : basis) {
    out << (*element.second);
  }
  out << std::flush;
  return out;
}

// adds a Gaussian of a pseudopotential
ECPGaussianPrimitive& ECPShell::addGaussian(Index power, double decay,
                                            double contraction) {
  _gaussians.push_back(ECPGaussianPrimitive(power, decay, contraction));
  return _gaussians.back();
}

}  // namespace xtp
}  // namespace votca
