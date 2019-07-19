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
#include "votca/xtp/ecpbasisset.h"
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

void ECPBasisSet::Load(const std::string& name) {
  tools::Property basis_property;
  _name = name;
  // if name contains .xml, assume a ecp .xml file is located in the working
  // directory
  std::size_t found_xml = name.find(".xml");
  std::string xmlFile;
  if (found_xml != std::string::npos) {
    xmlFile = name;
  } else {
    // get the path to the shared folders with xml files
    char* votca_share = getenv("VOTCASHARE");
    if (votca_share == NULL)
      throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/xtp/ecps/") +
              name + std::string(".xml");
  }
  bool success = load_property_from_xml(basis_property, xmlFile);

  if (!success) {
    throw std::runtime_error("ECP could not be loaded!");
  }

  std::vector<tools::Property*> elementProps =
      basis_property.Select("pseudopotential.element");

  for (tools::Property* elementProp : elementProps) {
    std::string elementName = elementProp->getAttribute<std::string>("name");
    int lmax = elementProp->getAttribute<int>("lmax");
    int ncore = elementProp->getAttribute<int>("ncore");

    ECPElement& element = addElement(elementName, lmax, ncore);

    std::vector<tools::Property*> shellProps = elementProp->Select("shell");
    for (tools::Property* shellProp : shellProps) {
      std::string shellType = shellProp->getAttribute<std::string>("type");
      if (shellType.size() > 1) {
        throw std::runtime_error(
            "In ecps no combined shells e.g. SP are allowed");
      }
      ECPShell& shell = element.addShell(shellType);
      std::vector<tools::Property*> constProps = shellProp->Select("constant");
      for (tools::Property* constProp : constProps) {
        int power = constProp->getAttribute<int>("power");
        double decay = constProp->getAttribute<double>("decay");
        double contraction = constProp->getAttribute<double>("contraction");
        shell.addGaussian(power, decay, contraction);
      }
    }
  }
  return;
}

// adding an Element to a Pseudopotential Library
ECPElement& ECPBasisSet::addElement(std::string elementType, int lmax,
                                    int ncore) {
  std::shared_ptr<ECPElement> element(new ECPElement(elementType, lmax, ncore));
  _elements[elementType] = element;
  return *element;
};

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

  out << "Type:" << shell.getType() << " Func: " << shell.getnumofFunc()
      << "\n";
  for (const auto& gaussian : shell._gaussians) {
    out << " Gaussian Decay: " << gaussian._decay;
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
ECPGaussianPrimitive& ECPShell::addGaussian(int power, double decay,
                                            double contraction) {
  _gaussians.push_back(ECPGaussianPrimitive(power, decay, contraction));
  return _gaussians.back();
}

}  // namespace xtp
}  // namespace votca
