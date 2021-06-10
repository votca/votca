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

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/tools/elements.h"

namespace votca {
namespace tools {

/*************************
 * Public Facing Methods *
 *************************/
bool Elements::isElement(std::string name) {
  return isEleShort(name) || isEleFull(name);
}

Index Elements::getNucCrg(std::string name) {
  if (!this->filled_NucCrg_) {
    this->FillNucCrg();
    filled_NucCrg_ = true;
  }
  try {
    return NucCrg_.at(name);
  } catch (const std::out_of_range& oor) {
    throw std::runtime_error("Nuclearcharge of element " + name +
                             " not found.");
  }
}

Index Elements::getEleNum(std::string name) {
  if (!this->filled_EleNum_) {
    this->FillEleNum();
    filled_EleNum_ = true;
  }
  try {
    return EleNum_.at(name);
  } catch (const std::out_of_range& oor) {
    throw std::runtime_error("Elementnumber of element " + name +
                             " not found.");
  }
}

double Elements::getMass(std::string name) {
  if (!this->filled_Mass_) {
    this->FillMass();
    filled_Mass_ = true;
  }

  try {
    return Mass_.at(name);
  } catch (const std::out_of_range& oor) {
    throw std::runtime_error("Mass of element " + name + " not found.");
  }
}

double Elements::getVdWChelpG(std::string name) {
  if (!this->filled_VdWChelpG_) {
    this->FillVdWChelpG();
    filled_VdWChelpG_ = true;
  }
  if (VdWChelpG_.count(name) == 0) {
    throw std::runtime_error("Element not found in VdWChelpG map " + name);
  }
  return VdWChelpG_.at(name);
}

double Elements::getVdWMK(std::string name) {
  if (!this->filled_VdWMK_) {
    this->FillVdWMK();
    filled_VdWMK_ = true;
  }
  if (VdWMK_.count(name) == 0) {
    throw std::runtime_error("Element not found in VdWMP map " + name);
  }
  return VdWMK_.at(name);
}

double Elements::getPolarizability(std::string name) {
  if (!this->filled_ElPolarizability_) {
    this->FillPolarizability();
    filled_ElPolarizability_ = true;
  }
  if (ElPolarizability_.count(name) == 0) {
    throw std::runtime_error("Element not found in ElPolarizability map " +
                             name);
  }
  return ElPolarizability_.at(name);
}

double Elements::getCovRad(std::string name, std::string unit) {
  // TODO - This should be replaced with an object, an object that should
  //       auto recognise the units and return it in a standard type
  if (!this->filled_CovRad_) {
    this->FillCovRad();
    filled_CovRad_ = true;
  }
  if (!unit.compare("bohr")) {
    return conv::ang2bohr * CovRad_.find(name)->second;
  }
  if (!unit.compare("nm")) {
    return conv::ang2nm * CovRad_.find(name)->second;
  }
  if (!unit.compare("ang")) {
    return CovRad_.find(name)->second;
  }
  throw std::runtime_error("Must specify appropriate units " + unit +
                           " is not known");
}

std::string Elements::getEleName(Index elenum) {
  if (!this->filled_EleName_) {
    this->FillEleName();
    filled_EleName_ = true;
  }
  return EleName_.at(elenum);
}

std::string Elements::getEleFull(std::string eleshort) {
  if (!this->filled_EleFull_) {
    this->FillEleFull();
    filled_EleFull_ = true;
  }
  return EleFull_.at(eleshort);
}

std::string Elements::getEleShort(std::string elefull) {
  if (!this->filled_EleShort_) {
    this->FillEleShort();
    filled_EleShort_ = true;
  }
  return EleShort_.at(elefull);
}

bool Elements::isEleFull(std::string fullname) {
  if (!this->filled_EleShort_) {
    this->FillEleShort();
    filled_EleShort_ = true;
  }
  std::string name_upper = boost::to_upper_copy<std::string>(fullname);
  return EleShort_.count(name_upper);
}

bool Elements::isEleShort(std::string shortname) {
  if (!this->filled_EleFull_) {
    this->FillEleFull();
    filled_EleFull_ = true;
  }
  return EleFull_.count(shortname);
}

bool Elements::isMassAssociatedWithElement(double mass, double tolerance) {
  auto closestMatch = findShortNameOfElementClosestInMass_(mass);
  if (closestMatch.second / Mass_[closestMatch.first] > tolerance) {
    return false;
  }
  return true;
}

std::string Elements::getEleShortClosestInMass(double mass, double tolerance) {
  auto closestMatch = findShortNameOfElementClosestInMass_(mass);
  if (closestMatch.second / Mass_[closestMatch.first] > tolerance) {
    throw std::runtime_error(
        "In attempt to determine if mass is associated "
        " with an element the mass exceeds tolerance of a possible match");
  }
  return closestMatch.first;
}

/*******************
 * Private Methods *
 *******************/

std::pair<std::string, double> Elements::findShortNameOfElementClosestInMass_(
    double mass) {
  if (!this->filled_Mass_) {
    this->FillMass();
    filled_Mass_ = true;
  }
  std::string eleShort = "H";
  double diff = std::fabs(mass - Mass_[eleShort]);
  for (const auto& ele_pr : Mass_) {
    if (std::fabs(ele_pr.second - mass) < diff) {
      eleShort = ele_pr.first;
      diff = std::fabs(ele_pr.second - mass);
    }
  }
  return std::pair<std::string, double>(eleShort, diff);
}

void Elements::FillMass() {
  // masses of atoms
  Mass_["H"] = 1.00794;
  Mass_["He"] = 4.002602;
  Mass_["Li"] = 6.941;
  Mass_["Be"] = 9.012182;
  Mass_["B"] = 10.811;
  Mass_["C"] = 12.0107;
  Mass_["N"] = 14.00674;
  Mass_["O"] = 15.9994;
  Mass_["F"] = 18.9984032;
  Mass_["Ne"] = 20.1797;
  Mass_["Na"] = 22.989770;
  Mass_["Mg"] = 24.3050;
  Mass_["Al"] = 26.981538;
  Mass_["Si"] = 28.0855;
  Mass_["P"] = 30.973761;
  Mass_["S"] = 32.066;
  Mass_["Cl"] = 35.4527;
  Mass_["Ar"] = 39.948;
  Mass_["K"] = 39.098;
  Mass_["Ca"] = 40.078;
  Mass_["Sc"] = 44.956;
  Mass_["Ti"] = 47.867;
  Mass_["V"] = 50.942;
  Mass_["Cr"] = 51.996;
  Mass_["Mn"] = 54.938;
  Mass_["Fe"] = 55.845;
  Mass_["Co"] = 58.933;
  Mass_["Ni"] = 58.693;
  Mass_["Cu"] = 63.546;
  Mass_["Zn"] = 65.38;
  Mass_["Ga"] = 69.723;
  Mass_["Ge"] = 72.630;
  Mass_["As"] = 74.922;
  Mass_["Se"] = 78.971;
  Mass_["Br"] = 79.90;
  Mass_["Kr"] = 83.798;
  Mass_["Rb"] = 85.468;
  Mass_["Sr"] = 87.62;
  Mass_["Y"] = 88.906;
  Mass_["Zr"] = 91.224;
  Mass_["Nb"] = 92.906;
  Mass_["Mo"] = 95.95;
  Mass_["Tc"] = 98.0;
  Mass_["Ru"] = 101.07;
  Mass_["Rh"] = 102.91;
  Mass_["Pd"] = 106.42;
  Mass_["Ag"] = 107.8682;
  Mass_["Cd"] = 112.41;
  Mass_["In"] = 114.82;
  Mass_["Sn"] = 118.71;
  Mass_["Sb"] = 121.76;
  Mass_["Te"] = 127.60;
  Mass_["I"] = 126.90;
  Mass_["Xe"] = 131.29;
  Mass_["Cs"] = 132.91;
  Mass_["Ba"] = 137.33;
  Mass_["Hf"] = 178.49;
  Mass_["Ta"] = 180.49;
  Mass_["W"] = 183.84;
  Mass_["Re"] = 186.21;
  Mass_["Os"] = 190.23;
  Mass_["Ir"] = 192.22;
  Mass_["Pt"] = 195.08;
  Mass_["Au"] = 196.97;
  Mass_["Hg"] = 200.59;
  Mass_["Tl"] = 204.38;
  Mass_["Pb"] = 207.2;
  Mass_["Bi"] = 208.98;
  Mass_["Po"] = 209;
  Mass_["At"] = 210;
  Mass_["Rn"] = 222;
}

void Elements::FillVdWChelpG() {
  // VdW radii in Angstrom as used in CHELPG paper [Journal of Computational
  // Chemistry 11, 361, 1990]and Gaussian
  VdWChelpG_["H"] = 1.45;
  VdWChelpG_["He"] = 1.45;
  VdWChelpG_["Li"] = 1.5;
  VdWChelpG_["Be"] = 1.5;
  VdWChelpG_["B"] = 1.5;
  VdWChelpG_["C"] = 1.5;
  VdWChelpG_["N"] = 1.7;
  VdWChelpG_["O"] = 1.7;
  VdWChelpG_["F"] = 1.7;
  VdWChelpG_["Ne"] = 1.7;
  VdWChelpG_["Na"] = 2.0;
  VdWChelpG_["Mg"] = 2.0;
  VdWChelpG_["Al"] = 2.0;
  VdWChelpG_["Si"] = 2.0;
  VdWChelpG_["P"] = 2.0;
  VdWChelpG_["S"] = 2.0;
  VdWChelpG_["Cl"] = 2.0;
  VdWChelpG_["Ar"] = 2.0;
  VdWChelpG_["Ag"] = 1.7;
}

void Elements::FillNucCrg() {
  // Nuclear Charges
  NucCrg_["H"] = 1.00;
  NucCrg_["He"] = 2.00;
  NucCrg_["Li"] = 3.00;
  NucCrg_["Be"] = 4.00;
  NucCrg_["B"] = 5.00;
  NucCrg_["C"] = 6.00;
  NucCrg_["N"] = 7.00;
  NucCrg_["O"] = 8.00;
  NucCrg_["F"] = 9.00;
  NucCrg_["Ne"] = 10.00;
  NucCrg_["Na"] = 11.00;
  NucCrg_["Mg"] = 12.00;
  NucCrg_["Al"] = 13.00;
  NucCrg_["Si"] = 14.00;
  NucCrg_["P"] = 15.00;
  NucCrg_["S"] = 16.00;
  NucCrg_["Cl"] = 17.00;
  NucCrg_["Ar"] = 18.00;
  NucCrg_["K"] = 19.00;
  NucCrg_["Ca"] = 20.00;
  NucCrg_["Sc"] = 21.00;
  NucCrg_["Ti"] = 22.00;
  NucCrg_["V"] = 23.00;
  NucCrg_["Cr"] = 24.00;
  NucCrg_["Mn"] = 25.00;
  NucCrg_["Fe"] = 26.00;
  NucCrg_["Co"] = 27.00;
  NucCrg_["Ni"] = 28.00;
  NucCrg_["Cu"] = 29.00;
  NucCrg_["Zn"] = 30.00;
  NucCrg_["Ga"] = 31.00;
  NucCrg_["Ge"] = 32.00;
  NucCrg_["As"] = 33.00;
  NucCrg_["Se"] = 34.00;
  NucCrg_["Br"] = 35.00;
  NucCrg_["Kr"] = 36.00;
  NucCrg_["Rb"] = 37.00;
  NucCrg_["Sr"] = 38.00;
  NucCrg_["Y"] = 39.00;
  NucCrg_["Zr"] = 40.00;
  NucCrg_["Nb"] = 41.00;
  NucCrg_["Mo"] = 42.00;
  NucCrg_["Tc"] = 43.00;
  NucCrg_["Ru"] = 44.00;
  NucCrg_["Rh"] = 45.00;
  NucCrg_["Pd"] = 46.00;
  NucCrg_["Ag"] = 47.00;
  NucCrg_["Cd"] = 48.00;
  NucCrg_["In"] = 49.00;
  NucCrg_["Sn"] = 50.00;
  NucCrg_["Sb"] = 51.00;
  NucCrg_["Te"] = 52.00;
  NucCrg_["I"] = 53.00;
  NucCrg_["Xe"] = 54.00;
  NucCrg_["Cs"] = 55.00;
  NucCrg_["Ba"] = 56.00;

  NucCrg_["Hf"] = 72.00;
  NucCrg_["Ta"] = 73.00;
  NucCrg_["W"] = 74.00;
  NucCrg_["Re"] = 75.00;
  NucCrg_["Os"] = 76.00;
  NucCrg_["Ir"] = 77.00;
  NucCrg_["Pt"] = 78.00;
  NucCrg_["Au"] = 79.00;
  NucCrg_["Hg"] = 80.00;
  NucCrg_["Tl"] = 81.00;
  NucCrg_["Pb"] = 82.00;
  NucCrg_["Bi"] = 83.00;
  NucCrg_["Po"] = 84.00;
  NucCrg_["At"] = 85.00;
  NucCrg_["Rn"] = 86.00;
}

void Elements::FillCovRad() {
  // Covalent Radii, used by BulkESP to break system into molecules
  // data from http://pubs.rsc.org/en/content/articlehtml/2008/dt/b801115j
  // values in [Angstroms]
  CovRad_["H"] = 0.31;
  CovRad_["He"] = 0.28;
  CovRad_["Li"] = 1.28;
  CovRad_["Be"] = 0.96;
  CovRad_["B"] = 0.84;
  // This is for sp3
  CovRad_["C"] = 0.76;
  CovRad_["N"] = 0.71;
  CovRad_["O"] = 0.66;
  CovRad_["F"] = 0.57;
  CovRad_["Ne"] = 0.58;
  CovRad_["Na"] = 1.66;
  CovRad_["Mg"] = 1.41;
  CovRad_["Al"] = 1.21;
  CovRad_["Si"] = 1.11;
  CovRad_["P"] = 1.07;
  CovRad_["S"] = 1.05;
  CovRad_["Cl"] = 1.02;
  CovRad_["Ar"] = 1.06;
  CovRad_["K"] = 2.03;
  CovRad_["Ca"] = 1.76;
  CovRad_["Sc"] = 1.70;
  CovRad_["Ti"] = 1.60;
  CovRad_["V"] = 1.53;
  CovRad_["Cr"] = 1.39;
  CovRad_["Mn"] = 1.61;
  CovRad_["Fe"] = 1.52;
  CovRad_["Co"] = 1.50;
  CovRad_["Ni"] = 1.24;
  CovRad_["Cu"] = 1.32;
  CovRad_["Zn"] = 1.22;
  CovRad_["Ga"] = 1.22;
  CovRad_["Ge"] = 1.20;
  CovRad_["As"] = 1.19;
  CovRad_["Se"] = 1.20;
  CovRad_["Br"] = 1.20;
  CovRad_["Kr"] = 1.16;
  CovRad_["Rb"] = 2.20;
  CovRad_["Sr"] = 1.95;
  CovRad_["Y"] = 1.90;
  CovRad_["Zr"] = 1.75;
  CovRad_["Nb"] = 1.64;
  CovRad_["Mo"] = 1.54;
  CovRad_["Tc"] = 1.47;
  CovRad_["Ru"] = 1.46;
  CovRad_["Rh"] = 1.42;
  CovRad_["Pd"] = 1.39;
  CovRad_["Ag"] = 1.45;
  CovRad_["Cd"] = 1.44;
  CovRad_["In"] = 1.42;
  CovRad_["Sn"] = 1.39;
  CovRad_["Sb"] = 1.39;
  CovRad_["Te"] = 1.38;
  CovRad_["I"] = 1.39;
  CovRad_["Xe"] = 1.40;
  CovRad_["Cs"] = 2.44;
  CovRad_["Ba"] = 2.15;
  CovRad_["Hf"] = 1.75;
  CovRad_["Ta"] = 1.70;
  CovRad_["W"] = 1.62;
  CovRad_["Re"] = 1.51;
  CovRad_["Os"] = 1.44;
  CovRad_["Ir"] = 1.41;
  CovRad_["Pt"] = 1.36;
  CovRad_["Au"] = 1.36;
  CovRad_["Hg"] = 1.32;
  CovRad_["Tl"] = 1.45;
  CovRad_["Pb"] = 1.46;
  CovRad_["Bi"] = 1.48;
  CovRad_["Po"] = 1.40;
  CovRad_["At"] = 1.50;
  CovRad_["Rn"] = 1.50;
}

void Elements::FillEleNum() {
  // Nuclear Charges
  EleNum_["H"] = 1;
  EleNum_["He"] = 2;
  EleNum_["Li"] = 3;
  EleNum_["Be"] = 4;
  EleNum_["B"] = 5;
  EleNum_["C"] = 6;
  EleNum_["N"] = 7;
  EleNum_["O"] = 8;
  EleNum_["F"] = 9;
  EleNum_["Ne"] = 10;
  EleNum_["Na"] = 11;
  EleNum_["Mg"] = 12;
  EleNum_["Al"] = 13;
  EleNum_["Si"] = 14;
  EleNum_["P"] = 15;
  EleNum_["S"] = 16;
  EleNum_["Cl"] = 17;
  EleNum_["Ar"] = 18;
  EleNum_["K"] = 19;
  EleNum_["Ca"] = 20;
  EleNum_["Sc"] = 21;
  EleNum_["Ti"] = 22;
  EleNum_["V"] = 23;
  EleNum_["Cr"] = 24;
  EleNum_["Mn"] = 25;
  EleNum_["Fe"] = 26;
  EleNum_["Co"] = 27;
  EleNum_["Ni"] = 28;
  EleNum_["Cu"] = 29;
  EleNum_["Zn"] = 30;
  EleNum_["Ga"] = 31;
  EleNum_["Ge"] = 32;
  EleNum_["As"] = 33;
  EleNum_["Se"] = 34;
  EleNum_["Br"] = 35;
  EleNum_["Kr"] = 36;
  EleNum_["Rb"] = 37;
  EleNum_["Sr"] = 38;
  EleNum_["Y"] = 39;
  EleNum_["Zr"] = 40;
  EleNum_["Nb"] = 41;
  EleNum_["Mo"] = 42;
  EleNum_["Tc"] = 43;
  EleNum_["Ru"] = 44;
  EleNum_["Rh"] = 45;
  EleNum_["Pd"] = 46;
  EleNum_["Ag"] = 47;
  EleNum_["Cd"] = 48;
  EleNum_["In"] = 49;
  EleNum_["Sn"] = 50;
  EleNum_["Sb"] = 51;
  EleNum_["Te"] = 52;
  EleNum_["I"] = 53;
  EleNum_["Xe"] = 54;
  EleNum_["Cs"] = 55;
  EleNum_["Ba"] = 56;
  EleNum_["Hf"] = 72;
  EleNum_["Ta"] = 73;
  EleNum_["W"] = 74;
  EleNum_["Re"] = 75;
  EleNum_["Os"] = 76;
  EleNum_["Ir"] = 77;
  EleNum_["Pt"] = 78;
  EleNum_["Au"] = 79;
  EleNum_["Hg"] = 80;
  EleNum_["Tl"] = 81;
  EleNum_["Pb"] = 82;
  EleNum_["Bi"] = 83;
  EleNum_["Po"] = 84;
  EleNum_["At"] = 85;
  EleNum_["Rn"] = 86;
}

void Elements::FillEleName() {
  // Nuclear Charges
  EleName_[1] = "H";
  EleName_[2] = "He";
  EleName_[3] = "Li";
  EleName_[4] = "Be";
  EleName_[5] = "B";
  EleName_[6] = "C";
  EleName_[7] = "N";
  EleName_[8] = "O";
  EleName_[9] = "F";
  EleName_[10] = "Ne";
  EleName_[11] = "Na";
  EleName_[12] = "Mg";
  EleName_[13] = "Al";
  EleName_[14] = "Si";
  EleName_[15] = "P";
  EleName_[16] = "S";
  EleName_[17] = "Cl";
  EleName_[18] = "Ar";
  EleName_[19] = "K";
  EleName_[20] = "Ca";
  EleName_[21] = "Sc";
  EleName_[22] = "Ti";
  EleName_[23] = "V";
  EleName_[24] = "Cr";
  EleName_[25] = "Mn";
  EleName_[26] = "Fe";
  EleName_[27] = "Co";
  EleName_[28] = "Ni";
  EleName_[29] = "Cu";
  EleName_[30] = "Zn";
  EleName_[31] = "Ga";
  EleName_[32] = "Ge";
  EleName_[33] = "As";
  EleName_[34] = "Se";
  EleName_[35] = "Br";
  EleName_[36] = "Kr";
  EleName_[37] = "Rb";
  EleName_[38] = "Sr";
  EleName_[39] = "Y";
  EleName_[40] = "Zr";
  EleName_[41] = "Nb";
  EleName_[42] = "Mo";
  EleName_[43] = "Tc";
  EleName_[44] = "Ru";
  EleName_[45] = "Rh";
  EleName_[46] = "Pd";
  EleName_[47] = "Ag";
  EleName_[48] = "Cd";
  EleName_[49] = "In";
  EleName_[50] = "Sn";
  EleName_[51] = "Sb";
  EleName_[52] = "Te";
  EleName_[53] = "I";
  EleName_[54] = "Xe";
  EleName_[55] = "Cs";
  EleName_[56] = "Ba";

  EleName_[72] = "Hf";
  EleName_[73] = "Ta";
  EleName_[74] = "W";
  EleName_[75] = "Re";
  EleName_[76] = "Os";
  EleName_[77] = "Ir";
  EleName_[78] = "Pt";
  EleName_[79] = "Au";
  EleName_[80] = "Hg";
  EleName_[81] = "Tl";
  EleName_[82] = "Pb";
  EleName_[83] = "Bi";
  EleName_[84] = "Po";
  EleName_[85] = "At";
  EleName_[86] = "Rn";
}

void Elements::FillEleShort() {
  // VdW radii in Angstrom as used in MK Gaussian
  EleShort_["HYDROGEN"] = "H";
  EleShort_["HELIUM"] = "He";
  EleShort_["LITHIUM"] = "Li";
  EleShort_["BERYLLIUM"] = "Be";
  EleShort_["BORON"] = "B";
  EleShort_["CARBON"] = "C";
  EleShort_["NITROGEN"] = "Ni";
  EleShort_["OXYGEN"] = "O";
  EleShort_["FLUORINE"] = "Fl";
  EleShort_["NEON"] = "Ne";
  EleShort_["SODIUM"] = "Na";
  EleShort_["MAGNESIUM"] = "Mg";
  EleShort_["ALUMINUM"] = "Al";
  EleShort_["SILICON"] = "Si";
  EleShort_["PHOSPHORUS"] = "Ph";
  EleShort_["SULFUR"] = "S";
  EleShort_["CLORINE"] = "Cl";
  EleShort_["ARGON"] = "Ar";
  EleShort_["POTASSIUM"] = "K";
  EleShort_["CALCIUM"] = "Ca";
  EleShort_["SCANDDIUM"] = "Sc";
  EleShort_["TITANIUM"] = "Ti";
  EleShort_["VANADIUM"] = "V";
  EleShort_["CHROMIUM"] = "Cr";
  EleShort_["MANGANESE"] = "Mn";
  EleShort_["IRON"] = "Fe";
  EleShort_["COBALT"] = "Co";
  EleShort_["NICKEL"] = "Ni";
  EleShort_["COPPER"] = "Cu";
  EleShort_["ZINC"] = "Zn";
  EleShort_["GALLIUM"] = "Ga";
  EleShort_["GERMANIUM"] = "Ge";
  EleShort_["ARSENIC"] = "As";
  EleShort_["SELENIUM"] = "Se";
  EleShort_["BROMINE"] = "Br";
  EleShort_["KRPTON"] = "Kr";
  EleShort_["RUBIDIUM"] = "Rb";
  EleShort_["STRONTIUM"] = "Sr";
  EleShort_["YTTRIUM"] = "Y";
  EleShort_["ZIRCONIUM"] = "Zr";
  EleShort_["NIOBIUM"] = "Nb";
  EleShort_["MOLYBDENUM"] = "Mo";
  EleShort_["TECHNETIUM"] = "Tc";
  EleShort_["RUTHENIUM"] = "Ru";
  EleShort_["RHODIUM"] = "Rh";
  EleShort_["PALLADIUM"] = "Pd";
  EleShort_["SILVER"] = "Ag";
  EleShort_["CADMIUM"] = "Cd";
  EleShort_["INDIUM"] = "In";
  EleShort_["TIN"] = "Sn";
  EleShort_["ANTIMONY"] = "Sb";
  EleShort_["TELLURIUM"] = "Te";
  EleShort_["IODINE"] = "I";
  EleShort_["XENON"] = "Xe";
  EleShort_["CAESIUM"] = "Cs";
  EleShort_["BARIUM"] = "Ba";
  EleShort_["HAFNIUM"] = "Hf";
  EleShort_["TANTALUM"] = "Ta";
  EleShort_["TUNGSTEN"] = "W";
  EleShort_["RHENIUM"] = "Re";
  EleShort_["OSMIUM"] = "Os";
  EleShort_["IRIDIUM"] = "Ir";
  EleShort_["PLATINUM"] = "Pt";
  EleShort_["GOLD"] = "Au";
  EleShort_["MERCURY"] = "Hg";
  EleShort_["THALLIUM"] = "Tl";
  EleShort_["LEAD"] = "Pd";
  EleShort_["BISMUTH"] = "Bi";
  EleShort_["PLONIUM"] = "Po";
  EleShort_["ASTATINE"] = "At";
  EleShort_["RADON"] = "Rn";
}

void Elements::FillEleFull() {
  // VdW radii in Angstrom as used in MK Gaussian
  EleFull_["H"] = "HYDROGEN";
  EleFull_["He"] = "HELIUM";
  EleFull_["Li"] = "LITHIUM";
  EleFull_["Be"] = "BERYLLIUM";
  EleFull_["B"] = "BORON";
  EleFull_["C"] = "CARBON";
  EleFull_["N"] = "NITROGEN";
  EleFull_["O"] = "OXYGEN";
  EleFull_["F"] = "FLUORINE";
  EleFull_["Ne"] = "NEON";
  EleFull_["Na"] = "SODIUM";
  EleFull_["Mg"] = "MAGNESIUM";
  EleFull_["Al"] = "ALUMINUM";
  EleFull_["Si"] = "SILICON";
  EleFull_["P"] = "PHOSPHORUS";
  EleFull_["S"] = "SULFUR";
  EleFull_["Cl"] = "CHLORINE";
  EleFull_["Ar"] = "ARGON";
  EleFull_["K"] = "POTASSIUM";
  EleFull_["Ca"] = "CALCIUM";
  EleFull_["Sc"] = "SCANDIUM";
  EleFull_["Ti"] = "TITANIUM";
  EleFull_["V"] = "VANADIUM";
  EleFull_["Cr"] = "CHROMIUM";
  EleFull_["Mn"] = "MANGANESE";
  EleFull_["Fe"] = "IRON";
  EleFull_["Co"] = "COBALT";
  EleFull_["Ni"] = "NICKEL";
  EleFull_["Cu"] = "COPPER";
  EleFull_["Zn"] = "ZINC";
  EleFull_["Ga"] = "GALLIUM";
  EleFull_["Ge"] = "GERMANIUM";
  EleFull_["As"] = "ARSENIC";
  EleFull_["Se"] = "SELENIUM";
  EleFull_["Br"] = "BROMINE";
  EleFull_["Kr"] = "KRYPTON";
  EleFull_["Rb"] = "RUBIDIUM";
  EleFull_["Sr"] = "STRONTIUM";
  EleFull_["Y"] = "YTTRIUM";
  EleFull_["Zr"] = "ZIRCONIUM";
  EleFull_["Nb"] = "NIOBIUM";
  EleFull_["Mo"] = "MOLYBDENUM";
  EleFull_["Tc"] = "TECHNETIUM";
  EleFull_["Ru"] = "TUTHENIUM";
  EleFull_["Rh"] = "RHODIUM";
  EleFull_["Pd"] = "PALLADIUM";
  EleFull_["Ag"] = "SILVER";
  EleFull_["Cd"] = "CADMIUM";
  EleFull_["In"] = "INDIUM";
  EleFull_["Sn"] = "TIN";
  EleFull_["Sb"] = "ANTIMONY";
  EleFull_["Te"] = "TELLURIUM";
  EleFull_["I"] = "IODINE";
  EleFull_["Xe"] = "XENON";
  EleFull_["Cs"] = "CEASIUM";
  EleFull_["Ba"] = "BARIUM";
  EleFull_["Hf"] = "HAFNIUM";
  EleFull_["Ta"] = "TANTALUM";
  EleFull_["W"] = "TUNGSTEN";
  EleFull_["Re"] = "RHENIUM";
  EleFull_["Os"] = "OSMIUM";
  EleFull_["Ir"] = "IRIDIUM";
  EleFull_["Pt"] = "PLATINUM";
  EleFull_["Au"] = "GOLD";
  EleFull_["Hg"] = "MERCURY";
  EleFull_["Tl"] = "THALLIUM";
  EleFull_["Pb"] = "LEAD";
  EleFull_["Bi"] = "BISMUTH";
  EleFull_["Po"] = "POLONIUM";
  EleFull_["At"] = "ASTATINE";
  EleFull_["Rn"] = "RADON";
}

void Elements::FillVdWMK() {
  // VdW radii in Angstrom as used in MK Gaussian
  VdWMK_["H"] = 1.2;
  VdWMK_["He"] = 1.2;
  VdWMK_["Li"] = 1.37;
  VdWMK_["Be"] = 1.45;
  VdWMK_["B"] = 1.45;
  VdWMK_["C"] = 1.5;
  VdWMK_["N"] = 1.5;
  VdWMK_["O"] = 1.4;
  VdWMK_["F"] = 1.35;
  VdWMK_["Ne"] = 1.3;
  VdWMK_["Na"] = 1.57;
  VdWMK_["Mg"] = 1.36;
  VdWMK_["Al"] = 1.24;
  VdWMK_["Si"] = 1.17;
  VdWMK_["P"] = 1.8;
  VdWMK_["S"] = 1.75;
  VdWMK_["Cl"] = 1.7;
  VdWMK_["Ag"] = 2.0;
}

void Elements::FillPolarizability() {
  ElPolarizability_["H"] = 0.496e-3;
  ElPolarizability_["C"] = 1.334e-3;
  ElPolarizability_["N"] = 1.073e-3;
  ElPolarizability_["O"] = 0.837e-3;
  ElPolarizability_["S"] = 2.926e-3;
  ElPolarizability_["F"] = 0.440e-3;
  ElPolarizability_["Si"] = 3.962e-3;  // B3LYP/6-311+g(2d,2p)
  ElPolarizability_["Zn"] = 5.962e-3;  // B3LYP/6-311+g(2d,2p)
  ElPolarizability_["Al"] =
      5.80e-3;  //[1]P. Fuentealba, “The static dipole polarizability of
                // aluminium atom: discrepancy between theory and experiment,”
                // Chemical physics letters, vol. 397, no. 4, pp. 459–461, 2004.
}

}  // namespace tools
}  // namespace votca
