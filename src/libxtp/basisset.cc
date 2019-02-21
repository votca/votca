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
#include "votca/xtp/basisset.h"
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

int FindLmax(const std::string& type) {
  int lmax;
  // single type shells
  if (type.length() == 1) {
    if (type == "S") {
      lmax = 0;
    } else if (type == "P") {
      lmax = 1;
    } else if (type == "D") {
      lmax = 2;
    } else if (type == "F") {
      lmax = 3;
    } else if (type == "G") {
      lmax = 4;
    } else if (type == "H") {
      lmax = 5;
    } else if (type == "I") {
      lmax = 6;
    } else {
      throw std::runtime_error("FindLmax: Shelltype not known");
    }
  } else {
    lmax = -1;
    for (unsigned i = 0; i < type.length(); ++i) {
      std::string local_shell = std::string(type, i, 1);
      int test = FindLmax(local_shell);
      if (test > lmax) {
        lmax = test;
      }
    }
  }
  return lmax;
}

int FindLmin(const std::string& type) {
  int lmin;
  if (type.length() == 1) {
    if (type == "S") {
      lmin = 0;
    } else if (type == "P") {
      lmin = 1;
    } else if (type == "D") {
      lmin = 2;
    } else if (type == "F") {
      lmin = 3;
    } else if (type == "G") {
      lmin = 4;
    } else if (type == "H") {
      lmin = 5;
    } else if (type == "I") {
      lmin = 6;
    }

    else {
      throw std::runtime_error("FindLmax: Shelltype not known");
    }
  } else {
    lmin = 1000;
    for (unsigned i = 0; i < type.length(); ++i) {
      std::string local_shell = std::string(type, i, 1);
      int test = FindLmin(local_shell);
      if (test == 0) {
        return 0;
      }
      if (test < lmin) {
        lmin = test;
      }
    }
  }
  return lmin;
}

int OffsetFuncShell(const std::string& shell_type) {
  int nbf;
  // single type shells
  if (shell_type.length() == 1) {
    if (shell_type == "S") {
      nbf = 0;
    } else if (shell_type == "P") {
      nbf = 1;
    } else if (shell_type == "D") {
      nbf = 4;
    } else if (shell_type == "F") {
      nbf = 9;
    } else if (shell_type == "G") {
      nbf = 16;
    } else if (shell_type == "H") {
      nbf = 25;
    } else if (shell_type == "I") {
      nbf = 36;
    } else {
      throw std::runtime_error("OffsetFuncShell: Shelltype not known");
    }
  } else {
    // for combined shells, go over all contributions and find minimal offset
    nbf = 1000;
    for (unsigned i = 0; i < shell_type.length(); ++i) {
      std::string local_shell = std::string(shell_type, i, 1);
      int test = OffsetFuncShell(local_shell);
      if (test < nbf) {
        nbf = test;
      }
    }
  }
  return nbf;
}

int NumFuncShell(const std::string& shell_type) {
  int nbf = 0;
  // single type shells
  if (shell_type.length() == 1) {
    if (shell_type == "S") {
      nbf = 1;
    } else if (shell_type == "P") {
      nbf = 3;
    } else if (shell_type == "D") {
      nbf = 5;
    } else if (shell_type == "F") {
      nbf = 7;
    } else if (shell_type == "G") {
      nbf = 9;
    } else if (shell_type == "H") {
      nbf = 11;
    } else if (shell_type == "I") {
      nbf = 13;
    } else {
      throw std::runtime_error("FindnumofFunc: Shelltype not known");
    }
  } else {
    // for combined shells, go over all contributions and add functions
    nbf = 0;
    for (unsigned i = 0; i < shell_type.length(); ++i) {
      std::string local_shell = std::string(shell_type, i, 1);
      nbf += NumFuncShell(local_shell);
    }
  }
  return nbf;
}

std::vector<int> NumFuncSubShell(const std::string& shell_type) {
  std::vector<int> subshells;
  // single type shells
  if (shell_type.length() == 1) {
    subshells.push_back(NumFuncShell(shell_type));
    // for combined shells, go over all contributions and add functions
  } else {
    for (unsigned i = 0; i < shell_type.length(); ++i) {
      std::string local_shell = std::string(shell_type, i, 1);
      subshells.push_back(NumFuncShell(local_shell));
    }
  }
  return subshells;
}

int NumFuncShell_cartesian(const std::string& shell_type) {
  int nbf;
  // single type shells defined here
  if (shell_type.length() == 1) {
    if (shell_type == "S") {
      nbf = 1;
    } else if (shell_type == "P") {
      nbf = 3;
    } else if (shell_type == "D") {
      nbf = 6;
    } else if (shell_type == "F") {
      nbf = 10;
    } else if (shell_type == "G") {
      nbf = 15;
    } else if (shell_type == "H") {
      nbf = 21;
    } else if (shell_type == "I") {
      nbf = 28;
    } else {
      throw std::runtime_error("NumFuncShell_cartesian shell_type not known");
    }
  } else {
    // for combined shells, sum over all contributions
    nbf = 0;
    for (unsigned i = 0; i < shell_type.length(); ++i) {
      std::string local_shell = std::string(shell_type, i, 1);
      nbf += NumFuncShell_cartesian(local_shell);
    }
  }

  return nbf;
}

int OffsetFuncShell_cartesian(const std::string& shell_type) {
  int nbf;
  // single type shells
  if (shell_type.length() == 1) {
    if (shell_type == "S") {
      nbf = 0;
    } else if (shell_type == "P") {
      nbf = 1;
    } else if (shell_type == "D") {
      nbf = 4;
    } else if (shell_type == "F") {
      nbf = 10;
    } else if (shell_type == "G") {
      nbf = 20;
    } else if (shell_type == "H") {
      nbf = 35;
    } else if (shell_type == "I") {
      nbf = 56;
    } else {
      throw std::runtime_error(
          "OffsetFuncShell_cartesian shell_type not known");
    }
  } else {
    // for combined shells, go over all contributions and find minimal offset
    nbf = 1000;
    for (unsigned i = 0; i < shell_type.length(); ++i) {
      std::string local_shell = std::string(shell_type, i, 1);
      int test = OffsetFuncShell_cartesian(local_shell);
      if (test < nbf) {
        nbf = test;
      }
    }
  }
  return nbf;
}

void BasisSet::LoadBasisSet(const std::string& name) {
  tools::Property basis_property;
  _name = name;
  // if name contains .xml, assume a basisset .xml file is located in the
  // working directory
  std::size_t found_xml = name.find(".xml");
  std::string xmlFile;
  if (found_xml != std::string::npos) {
    xmlFile = name;
  } else {
    // get the path to the shared folders with xml files
    char* votca_share = getenv("VOTCASHARE");
    if (votca_share == NULL)
      throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    xmlFile = std::string(getenv("VOTCASHARE")) +
              std::string("/xtp/basis_sets/") + name + std::string(".xml");
  }
  bool success = load_property_from_xml(basis_property, xmlFile);

  if (!success) {
    ;
  }
  std::vector<tools::Property*> elementProps =
      basis_property.Select("basis.element");

  for (tools::Property* elementProp : elementProps) {
    std::string elementName = elementProp->getAttribute<std::string>("name");
    Element& element = addElement(elementName);
    std::vector<tools::Property*> shellProps = elementProp->Select("shell");
    for (tools::Property* shellProp : shellProps) {
      std::string shellType = shellProp->getAttribute<std::string>("type");
      double shellScale = shellProp->getAttribute<double>("scale");

      Shell& shell = element.addShell(shellType, shellScale);
      std::vector<tools::Property*> constProps = shellProp->Select("constant");
      for (tools::Property* constProp : constProps) {
        double decay = constProp->getAttribute<double>("decay");
        std::vector<double> contraction =
            std::vector<double>(shell.getLmax() + 1, 0.0);
        std::vector<tools::Property*> contrProps =
            constProp->Select("contractions");
        for (tools::Property* contrProp : contrProps) {
          std::string contrType = contrProp->getAttribute<std::string>("type");
          double contrFactor = contrProp->getAttribute<double>("factor");
          if (contrType == "S")
            contraction[0] = contrFactor;
          else if (contrType == "P")
            contraction[1] = contrFactor;
          else if (contrType == "D")
            contraction[2] = contrFactor;
          else if (contrType == "F")
            contraction[3] = contrFactor;
          else if (contrType == "G")
            contraction[4] = contrFactor;
          else if (contrType == "H")
            contraction[5] = contrFactor;
          else if (contrType == "I")
            contraction[6] = contrFactor;
          else {
            throw std::runtime_error("LoadBasiset:Contractiontype not known");
          }
        }
        shell.addGaussian(decay, contraction);
      }
    }
  }
  return;
}

void BasisSet::LoadPseudopotentialSet(const std::string& name) {
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
    throw std::runtime_error("Basisset could not be loaded!");
  }

  std::vector<tools::Property*> elementProps =
      basis_property.Select("pseudopotential.element");

  for (tools::Property* elementProp : elementProps) {
    std::string elementName = elementProp->getAttribute<std::string>("name");
    int lmax = elementProp->getAttribute<int>("lmax");
    int ncore = elementProp->getAttribute<int>("ncore");

    Element& element = addElement(elementName, lmax, ncore);

    std::vector<tools::Property*> shellProps = elementProp->Select("shell");
    for (tools::Property* shellProp : shellProps) {
      std::string shellType = shellProp->getAttribute<std::string>("type");
      double shellScale = 1.0;

      Shell& shell = element.addShell(shellType, shellScale);

      std::vector<tools::Property*> constProps = shellProp->Select("constant");
      for (tools::Property* constProp : constProps) {
        int power = constProp->getAttribute<int>("power");
        double decay = constProp->getAttribute<double>("decay");
        std::vector<double> contraction;
        contraction.push_back(constProp->getAttribute<double>("contraction"));
        shell.addGaussian(power, decay, contraction);
      }
    }
  }
  return;
}

// adding an Element to a Basis Set
Element& BasisSet::addElement(std::string elementType) {
  std::shared_ptr<Element> element(new Element(elementType));
  _elements[elementType] = element;
  return *element;
};

// adding an Element to a Pseudopotential Library
Element& BasisSet::addElement(std::string elementType, int lmax, int ncore) {
  std::shared_ptr<Element> element(new Element(elementType, lmax, ncore));
  _elements[elementType] = element;
  return *element;
};

const Element& BasisSet::getElement(std::string element_type) const {
  std::map<std::string, std::shared_ptr<Element> >::const_iterator itm =
      _elements.find(element_type);
  if (itm == _elements.end()) {
    throw std::runtime_error("Basis set " + _name +
                             " does not have element of type " + element_type);
  }
  const Element& element = *((*itm).second);
  return element;
}

std::ostream& operator<<(std::ostream& out, const Shell& shell) {

  out << "Type:" << shell.getType() << " Scale:" << shell.getScale()
      << " Func: " << shell.getnumofFunc() << "\n";
  for (const auto& gaussian : shell._gaussians) {
    out << " Gaussian Decay: " << gaussian._decay;
    out << " Contractions:";
    for (const double& contraction : gaussian._contraction) {
      out << " " << contraction;
    }
    out << "\n";
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const Element& element) {
  out << "Element:" << element.getType() << "\n";
  for (const auto& shell : element) {
    out << shell;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const BasisSet& basis) {
  out << "BasisSet:" << basis._name << "\n";
  for (const auto& element : basis) {
    out << (*element.second);
  }
  out << std::flush;
  return out;
}

GaussianPrimitive& Shell::addGaussian(double decay,
                                      std::vector<double> contraction) {
  _gaussians.push_back(GaussianPrimitive(decay, contraction));
  return _gaussians.back();
}

// adds a Gaussian of a pseudopotential
GaussianPrimitive& Shell::addGaussian(int power, double decay,
                                      std::vector<double> contraction) {
  _gaussians.push_back(GaussianPrimitive(power, decay, contraction));
  return _gaussians.back();
}

}  // namespace xtp
}  // namespace votca
