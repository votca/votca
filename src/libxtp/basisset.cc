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
#include "votca/xtp/basisset.h"

namespace votca {
namespace xtp {

L StringToEnum(const std::string& type) {
  assert(!type.empty() && "Shelltype must be non empty!");
  assert(type.size() == 1 && "Shelltype size must be one");
  const char t = type.back();
  return StringToEnum(t);
}

L StringToEnum(char type) {
  L l;
  if (type == 'S') {
    l = L::S;
  } else if (type == 'P') {
    l = L::P;
  } else if (type == 'D') {
    l = L::D;
  } else if (type == 'F') {
    l = L::F;
  } else if (type == 'G') {
    l = L::G;
  } else if (type == 'H') {
    l = L::H;
  } else if (type == 'I') {
    l = L::I;
  } else {
    throw std::runtime_error("FindLmax: Shelltype '" + std::string(1, type) +
                             "' not known");
  }
  return l;
}

std::string EnumToString(L l) {
  switch (l) {
    case L::S:
      return "S";
    case L::P:
      return "P";
    case L::D:
      return "D";
    case L::F:
      return "F";
    case L::G:
      return "G";
    case L::H:
      return "H";
    case L::I:
      return "I";
  }
  return "";
}

Index OffsetFuncShell(L l) {
  switch (l) {
    case L::S:
      return 0;
    case L::P:
      return 1;
    case L::D:
      return 4;
    case L::F:
      return 9;
    case L::G:
      return 16;
    case L::H:
      return 25;
    case L::I:
      return 36;
  }
  return -1;
}

Index NumFuncShell(L l) { return 2 * Index(l) + 1; }

Index NumFuncShell_cartesian(L l) {
  Index lindex = Index(l);
  return (lindex + 1) * (lindex + 2) / 2;
}

Index OffsetFuncShell_cartesian(L l) {
  switch (l) {
    case L::S:
      return 0;
    case L::P:
      return 1;
    case L::D:
      return 4;
    case L::F:
      return 10;
    case L::G:
      return 20;
    case L::H:
      return 35;
    case L::I:
      return 56;
  }
  return -1;
}

bool CheckShellType(const std::string& shelltype) {
  if (shelltype.empty()) {
    return false;
  }
  std::vector<char> allowed_shells = {'S', 'P', 'D', 'F', 'G', 'H', 'I'};
  std::vector<char>::iterator it =
      std::find(allowed_shells.begin(), allowed_shells.end(), shelltype[0]);
  if (it == allowed_shells.end()) {
    return false;
  } else {
    Index index = std::distance(allowed_shells.begin(), it);
    for (Index i = 1; i < Index(shelltype.size()); i++) {
      if (index + i > Index(allowed_shells.size()) ||
          shelltype[i] != allowed_shells[index + i]) {
        return false;
      }
    }
  }

  return true;
}

void BasisSet::Load(const std::string& name) {

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
    if (votca_share == nullptr) {
      throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    }
    xmlFile = std::string(getenv("VOTCASHARE")) +
              std::string("/xtp/basis_sets/") + name + std::string(".xml");
  }
  tools::Property basis_property;
  basis_property.LoadFromXML(xmlFile);
  std::vector<tools::Property*> elementProps =
      basis_property.Select("basis.element");

  for (tools::Property* elementProp : elementProps) {
    std::string elementName = elementProp->getAttribute<std::string>("name");
    Element& element = addElement(elementName);
    std::vector<tools::Property*> shellProps = elementProp->Select("shell");
    for (tools::Property* shellProp : shellProps) {
      std::string shellType = shellProp->getAttribute<std::string>("type");
      if (!CheckShellType(shellType)) {
        throw std::runtime_error("Shelltype: '" + shellType +
                                 "' is not a valid shelltype!");
      }
      for (char subtype : shellType) {

        double shellScale = shellProp->getAttribute<double>("scale");

        Shell& shell = element.addShell(StringToEnum(subtype), shellScale);
        std::vector<tools::Property*> constProps =
            shellProp->Select("constant");
        for (tools::Property* constProp : constProps) {
          double decay = constProp->getAttribute<double>("decay");
          std::vector<tools::Property*> contrProps =
              constProp->Select("contractions");
          double contraction = 0.0;
          for (tools::Property* contrProp : contrProps) {
            std::string contrType =
                contrProp->getAttribute<std::string>("type");
            if (contrType != std::string(1, subtype)) {
              continue;
            }
            contraction = contrProp->getAttribute<double>("factor");
          }
          shell.addGaussian(decay, contraction);
        }
      }
    }
  }
  return;
}

// adding an Element to a Basis Set
Element& BasisSet::addElement(std::string elementType) {
  auto e = _elements.insert({elementType, Element(elementType)});
  if (!e.second) {
    throw std::runtime_error("Inserting element into basisset failed!");
  }
  return e.first->second;
}

const Element& BasisSet::getElement(std::string element_type) const {
  std::map<std::string, Element>::const_iterator itm =
      _elements.find(element_type);
  if (itm == _elements.end()) {
    throw std::runtime_error("Basis set " + _name +
                             " does not have element of type " + element_type);
  }
  return itm->second;
}

std::ostream& operator<<(std::ostream& out, const Shell& shell) {

  out << "Type:" << EnumToString(shell.getL()) << " Scale:" << shell.getScale()
      << " Func: " << shell.getnumofFunc() << "\n";
  for (const auto& gaussian : shell._gaussians) {
    out << " Gaussian Decay: " << gaussian.decay();
    out << " Contraction: " << gaussian.contraction();
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
    out << element.second;
  }
  out << std::flush;
  return out;
}

GaussianPrimitive& Shell::addGaussian(double decay, double contraction) {
  _gaussians.push_back(GaussianPrimitive(decay, contraction));
  return _gaussians.back();
}

}  // namespace xtp
}  // namespace votca
