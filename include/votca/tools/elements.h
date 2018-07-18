/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __VOTCA_TOOLS_ELEMENTS_H
#define __VOTCA_TOOLS_ELEMENTS_H

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <votca/tools/constants.h>

namespace votca {
namespace tools {

/**
    \brief information about an element

    The Elements class stores properties of several elements in the periodic
    table includeing: mass, element name, element symbol, van der Waals radius,
    effective nuclear potential, element number as it appears in the periodic
    table.

*/
class Elements {
 public:
  Elements() { FillMaps(); };
  ~Elements(){};

  /// ChelpG is a method for calculating the electrostatic potential outside
  /// of a molecule. The VdWChelpG radius in this case is specific to each
  /// element and defines the radius outside of which the electrostatic
  /// potential is calculated. The electrostatic potential within the radius
  /// is not calculated. For more information see the reference paper CHELPG
  /// paper [Journal of Computational Chemistry 11, 361, 1990]
  const double &getVdWChelpG(std::string name) const {
    if (_VdWChelpG.count(name) == 0)
      throw std::invalid_argument("Element not found in VdWChelpG map " + name);
    return _VdWChelpG.at(name);
  }

  /// Merz-Singh-Kollman MK method is a similar method for calculating the
  /// electrostatic potential outside of a molecule. Details of the method can
  /// be found here [Singh, U. Chandra (04/01/1984). "An approach to computing
  /// electrostatic charges for molecules". Journal of computational chemistry
  /// 5 (2), 129, 1984]. The VdWMK method will return the radii used to define
  /// where the electrostatic grid starts.
  const double &getVdWMK(std::string name) const {
    if (_VdWMK.count(name) == 0)
      throw std::invalid_argument("Element not found in VdWMP map " + name);
    return _VdWMK.at(name);
  }

  /// This method is now deprecated but was originally used for determining
  /// the effective nuclear potential of an atom. There are now more
  /// appropriate methods for doing this in the basisset.h and aobasis.h file
  /// in xtp.
  const double &getNucCrgECP(std::string name) const {
    throw std::invalid_argument("CrgECP map is deprecated");
  }

  /// Return the Nuclear charges of each atom. H - 1, He - 2, Na - 3 etc...
  const double &getNucCrg(std::string name) const { return _NucCrg.at(name); }
  /// Similar to the Nuclear charges but returns in integer form represents
  /// the id of the atom in the periodic table, the id starts at 1
  const int &getEleNum(std::string name) const { return _EleNum.at(name); }
  /// Returns the mass of each atom in a.u.
  const double &getMass(std::string name) const { return _Mass.at(name); }
  /// Returns the atomic polarisability of atom 
  // All polarizabilities in nm**3
  // Isotropic polarizability volume is evaluated from the tensor
  // as (a_xx * a_yy * a_zz )^(1/3) for eigenvalues of the polarizability tensor
  const double &getPolarizability(std::string name) const { if (_ElPolarizability.count(name) == 0)
      throw std::invalid_argument("Element not found in ElPolarizability map " + name);
    return _ElPolarizability.at(name); }
  /// Returns the covalent Radii of the atom 
  double getCovRad(std::string name,std::string unit ) const 
  { 
    //TODO - This should be replaced with an object, an object that should
    //       auto recognise the units and return it in a standard type
    if(!unit.compare("bohr")) return votca::tools::conv::ang2bohr*_CovRad.find(name)->second;
    if(!unit.compare("nm"))   return votca::tools::conv::ang2nm*_CovRad.find(name)->second;
    if(!unit.compare("ang"))  return _CovRad.find(name)->second; 
    
    throw std::invalid_argument("Must specify appropriate units " + unit + " is not known");
  }
  /// Provided the element number returns the symbol for the element name
  /// (1) = "H", (2) = "He", ...
  const std::string &getEleName(int elenum) const {
    return _EleName.at(elenum);
  }
  /// Provided the full element name returns the element symbol
  /// "Hydrogen" = "H", "HELIUM" = "He",...
  const std::string &getEleShort(std::string elefull) const {
    return _EleShort.at(elefull);
  }

	/// Get the shortened element name given a mass similar in size to one of 
	/// the elements. 
	std::string getEleShortClosestInMass(double mass){
		std::string eleShort = "H";
		double diff = abs(mass-_Mass[eleShort]);
		for(auto & ele_pr : _Mass ){
			if(abs(ele_pr.second-mass)<diff){
				eleShort = ele_pr.first;
				diff = abs(ele_pr.second-mass);
			}
		}
		return eleShort;
	}

  /// Provided the element symbol returns the element name
  /// "Pb" = "LEAD", "Na" = "SODIUM", ....
  const std::string &getEleFull(std::string eleshort) const {
    return _EleFull.at(eleshort);
  }

 private:
  std::map<std::string, double> _VdWChelpG;
  std::map<std::string, double> _VdWMK;
  std::map<std::string, double> _NucCrgECP;
  std::map<std::string, double> _NucCrg;
  std::map<std::string, double> _CovRad;
  std::map<std::string, double> _Mass;
  std::map<std::string, int>    _EleNum;
  std::map<int, std::string>    _EleName;
  
  std::map<std::string,double> _ElPolarizability;

  std::map<std::string, std::string> _EleShort;
  std::map<std::string, std::string> _EleFull;

  inline void FillMaps() {

    FillVdWChelpG();
    FillVdWMK();
    FillCovRad();
    FillNucCrgECP();
    FillNucCrg();
    FillEleNum();
    FillEleName();
    FillEleShort();
    FillEleFull();
    FillMass();
    FillPolarizability();
  };

  inline void FillMass() {
    // masses of atoms
    _Mass["H"]  = 1.00794;
    _Mass["He"] = 4.002602;
    _Mass["Li"] = 6.941;
    _Mass["Be"] = 9.012182;
    _Mass["B"]  = 10.811;
    _Mass["C"]  = 12.0107;
    _Mass["N"]  = 14.00674;
    _Mass["O"]  = 15.9994;
    _Mass["F"]  = 18.9984032;
    _Mass["Ne"] = 20.1797;
    _Mass["Na"] = 22.989770;
    _Mass["Mg"] = 24.3050;
    _Mass["Al"] = 26.981538;
    _Mass["Si"] = 28.0855;
    _Mass["P"]  = 30.973761;
    _Mass["S"]  = 32.066;
    _Mass["Cl"] = 35.4527;
    _Mass["Ar"] = 39.948;
    _Mass["K"]  = 39.098;
    _Mass["Ca"] = 40.078;
    _Mass["Sc"] = 44.956;
    _Mass["Ti"] = 47.867;
    _Mass["V"]  = 50.942;
    _Mass["Cr"] = 51.996;
    _Mass["Mn"] = 54.938;
    _Mass["Fe"] = 55.845;
    _Mass["Co"] = 58.933;
    _Mass["Ni"] = 58.693;
    _Mass["Cu"] = 63.546;
    _Mass["Zn"] = 65.38;
    _Mass["Ga"] = 69.723;
    _Mass["Ge"] = 72.630;
    _Mass["As"] = 74.922;
    _Mass["Se"] = 78.971;
    _Mass["Br"] = 79.90;
    _Mass["Kr"] = 83.798;
    _Mass["Rb"] = 85.468;
    _Mass["Sr"] = 87.62;
    _Mass["Y"]  = 88.906;
    _Mass["Zr"] = 91.224;
    _Mass["Nb"] = 92.906;
    _Mass["Mo"] = 95.95;
    _Mass["Tc"] = 98.0;
    _Mass["Ru"] = 101.07;
    _Mass["Rh"] = 102.91;
    _Mass["Pd"] = 106.42;
    _Mass["Ag"] = 107.8682;
    _Mass["Cd"] = 112.41;
    _Mass["In"] = 114.82;
    _Mass["Sn"] = 118.71;
    _Mass["Sb"] = 121.76;
    _Mass["Te"] = 127.60;
    _Mass["I"]  = 126.90;
    _Mass["Xe"] = 131.29;
    _Mass["Cs"] = 132.91;
    _Mass["Ba"] = 137.33;
    _Mass["Hf"] = 178.49;
    _Mass["Ta"] = 180.49;
    _Mass["W"]  = 183.84;
    _Mass["Re"] = 186.21;
    _Mass["Os"] = 190.23;
    _Mass["Ir"] = 192.22;
    _Mass["Pt"] = 195.08;
    _Mass["Au"] = 196.97;
    _Mass["Hg"] = 200.59;
    _Mass["Tl"] = 204.38;
    _Mass["Pb"] = 207.2;
    _Mass["Bi"] = 208.98;
    _Mass["Po"] = 209;
    _Mass["At"] = 210;
    _Mass["Rn"] = 222;
  };

  inline void FillVdWChelpG() {
    // VdW radii in Angstrom as used in CHELPG paper [Journal of Computational
    // Chemistry 11, 361, 1990]and Gaussian
    _VdWChelpG["H"]  = 1.45;
    _VdWChelpG["He"] = 1.45;
    _VdWChelpG["Li"] = 1.5;
    _VdWChelpG["Be"] = 1.5;
    _VdWChelpG["B"]  = 1.5;
    _VdWChelpG["C"]  = 1.5;
    _VdWChelpG["N"]  = 1.7;
    _VdWChelpG["O"]  = 1.7;
    _VdWChelpG["F"]  = 1.7;
    _VdWChelpG["Ne"] = 1.7;
    _VdWChelpG["Na"] = 2.0;
    _VdWChelpG["Mg"] = 2.0;
    _VdWChelpG["Al"] = 2.0;
    _VdWChelpG["Si"] = 2.0;
    _VdWChelpG["P"]  = 2.0;
    _VdWChelpG["S"]  = 2.0;
    _VdWChelpG["Cl"] = 2.0;
    _VdWChelpG["Ar"] = 2.0;
    _VdWChelpG["Ag"] = 1.7;
  };

  inline void FillNucCrgECP() {
    // Effective Core Potential Nuclear Charges
    _NucCrgECP["H"]  = 1.00;
    _NucCrgECP["He"] = 2.00;
    _NucCrgECP["Li"] = 1.00;
    _NucCrgECP["Be"] = 2.00;
    _NucCrgECP["B"]  = 3.00;
    _NucCrgECP["C"]  = 4.00;
    _NucCrgECP["N"]  = 5.00;
    _NucCrgECP["O"]  = 6.00;
    _NucCrgECP["F"]  = 7.00;
    _NucCrgECP["Ne"] = 8.00;
    _NucCrgECP["Na"] = 1.00;
    _NucCrgECP["Mg"] = 2.00;
    _NucCrgECP["Al"] = 3.00;
    _NucCrgECP["Si"] = 4.00;
    _NucCrgECP["P"]  = 5.00;
    _NucCrgECP["S"]  = 6.00;
    _NucCrgECP["Cl"] = 7.00;
    _NucCrgECP["Ar"] = 8.00;
    _NucCrgECP["Ag"] = 19.00;
  };

  inline void FillNucCrg() {
    // Nuclear Charges
    _NucCrg["H"]  = 1.00;
    _NucCrg["He"] = 2.00;
    _NucCrg["Li"] = 3.00;
    _NucCrg["Be"] = 4.00;
    _NucCrg["B"]  = 5.00;
    _NucCrg["C"]  = 6.00;
    _NucCrg["N"]  = 7.00;
    _NucCrg["O"]  = 8.00;
    _NucCrg["F"]  = 9.00;
    _NucCrg["Ne"] = 10.00;
    _NucCrg["Na"] = 11.00;
    _NucCrg["Mg"] = 12.00;
    _NucCrg["Al"] = 13.00;
    _NucCrg["Si"] = 14.00;
    _NucCrg["P"]  = 15.00;
    _NucCrg["S"]  = 16.00;
    _NucCrg["Cl"] = 17.00;
    _NucCrg["Ar"] = 18.00;
    _NucCrg["K"]  = 19.00;
    _NucCrg["Ca"] = 20.00;
    _NucCrg["Sc"] = 21.00;
    _NucCrg["Ti"] = 22.00;
    _NucCrg["V"]  = 23.00;
    _NucCrg["Cr"] = 24.00;
    _NucCrg["Mn"] = 25.00;
    _NucCrg["Fe"] = 26.00;
    _NucCrg["Co"] = 27.00;
    _NucCrg["Ni"] = 28.00;
    _NucCrg["Cu"] = 29.00;
    _NucCrg["Zn"] = 30.00;
    _NucCrg["Ga"] = 31.00;
    _NucCrg["Ge"] = 32.00;
    _NucCrg["As"] = 33.00;
    _NucCrg["Se"] = 34.00;
    _NucCrg["Br"] = 35.00;
    _NucCrg["Kr"] = 36.00;
    _NucCrg["Rb"] = 37.00;
    _NucCrg["Sr"] = 38.00;
    _NucCrg["Y"]  = 39.00;
    _NucCrg["Zr"] = 40.00;
    _NucCrg["Nb"] = 41.00;
    _NucCrg["Mo"] = 42.00;
    _NucCrg["Tc"] = 43.00;
    _NucCrg["Ru"] = 44.00;
    _NucCrg["Rh"] = 45.00;
    _NucCrg["Pd"] = 46.00;
    _NucCrg["Ag"] = 47.00;
    _NucCrg["Cd"] = 48.00;
    _NucCrg["In"] = 49.00;
    _NucCrg["Sn"] = 50.00;
    _NucCrg["Sb"] = 51.00;
    _NucCrg["Te"] = 52.00;
    _NucCrg["I"]  = 53.00;
    _NucCrg["Xe"] = 54.00;
    _NucCrg["Cs"] = 55.00;
    _NucCrg["Ba"] = 56.00;

    _NucCrg["Hf"] = 72.00;
    _NucCrg["Ta"] = 73.00;
    _NucCrg["W"]  = 74.00;
    _NucCrg["Re"] = 75.00;
    _NucCrg["Os"] = 76.00;
    _NucCrg["Ir"] = 77.00;
    _NucCrg["Pt"] = 78.00;
    _NucCrg["Au"] = 79.00;
    _NucCrg["Hg"] = 80.00;
    _NucCrg["Tl"] = 81.00;
    _NucCrg["Pb"] = 82.00;
    _NucCrg["Bi"] = 83.00;
    _NucCrg["Po"] = 84.00;
    _NucCrg["At"] = 85.00;
    _NucCrg["Rn"] = 86.00;
  };


  inline void FillCovRad(){
  // Covalent Radii, used by BulkESP to break system into molecules
  //data from http://pubs.rsc.org/en/content/articlehtml/2008/dt/b801115j
  //values in [Angstroms]
    _CovRad["H"]  = 0.31; 
    _CovRad["He"] = 0.28;
    _CovRad["Li"] = 1.28;
    _CovRad["Be"] = 0.96;
    _CovRad["B"]  = 0.84;
  // This is for sp3 
    _CovRad["C"]  = 0.76;
    _CovRad["N"]  = 0.71;
    _CovRad["O"]  = 0.66;
    _CovRad["F"]  = 0.57;
    _CovRad["Ne"] = 0.58;
    _CovRad["Na"] = 1.66;
    _CovRad["Mg"] = 1.41;
    _CovRad["Al"] = 1.21;
    _CovRad["Si"] = 1.11;
    _CovRad["P"]  = 1.07;
    _CovRad["S"]  = 1.05;
    _CovRad["Cl"] = 1.02;
    _CovRad["Ar"] = 1.06;
    _CovRad["K"]  = 2.03;
    _CovRad["Ca"] = 1.76;
    _CovRad["Sc"] = 1.70;
    _CovRad["Ti"] = 1.60;
    _CovRad["V"]  = 1.53;
    _CovRad["Cr"] = 1.39;
    _CovRad["Mn"] = 1.61;
    _CovRad["Fe"] = 1.52;
    _CovRad["Co"] = 1.50;
    _CovRad["Ni"] = 1.24;
    _CovRad["Cu"] = 1.32;
    _CovRad["Zn"] = 1.22;
    _CovRad["Ga"] = 1.22;
    _CovRad["Ge"] = 1.20;
    _CovRad["As"] = 1.19;
    _CovRad["Se"] = 1.20;
    _CovRad["Br"] = 1.20;
    _CovRad["Kr"] = 1.16;
    _CovRad["Rb"] = 2.20;
    _CovRad["Sr"] = 1.95;
    _CovRad["Y"]  = 1.90;
    _CovRad["Zr"] = 1.75;
    _CovRad["Nb"] = 1.64;
    _CovRad["Mo"] = 1.54;
    _CovRad["Tc"] = 1.47;
    _CovRad["Ru"] = 1.46;
    _CovRad["Rh"] = 1.42;
    _CovRad["Pd"] = 1.39;
    _CovRad["Ag"] = 1.45;
    _CovRad["Cd"] = 1.44;
    _CovRad["In"] = 1.42;
    _CovRad["Sn"] = 1.39;
    _CovRad["Sb"] = 1.39;
    _CovRad["Te"] = 1.38;
    _CovRad["I"]  = 1.39;
    _CovRad["Xe"] = 1.40;
    _CovRad["Cs"] = 2.44;
    _CovRad["Ba"] = 2.15;
    _CovRad["Hf"] = 1.75;
    _CovRad["Ta"] = 1.70;
    _CovRad["W"]  = 1.62;
    _CovRad["Re"] = 1.51;
    _CovRad["Os"] = 1.44;
    _CovRad["Ir"] = 1.41;
    _CovRad["Pt"] = 1.36;
    _CovRad["Au"] = 1.36;
    _CovRad["Hg"] = 1.32;
    _CovRad["Tl"] = 1.45;
    _CovRad["Pb"] = 1.46;
    _CovRad["Bi"] = 1.48;
    _CovRad["Po"] = 1.40;
    _CovRad["At"] = 1.50;
    _CovRad["Rn"] = 1.50;

  };

  inline void FillEleNum() {
    // Nuclear Charges
    _EleNum["H"]  = 1;
    _EleNum["He"] = 2;
    _EleNum["Li"] = 3;
    _EleNum["Be"] = 4;
    _EleNum["B"]  = 5;
    _EleNum["C"]  = 6;
    _EleNum["N"]  = 7;
    _EleNum["O"]  = 8;
    _EleNum["F"]  = 9;
    _EleNum["Ne"] = 10;
    _EleNum["Na"] = 11;
    _EleNum["Mg"] = 12;
    _EleNum["Al"] = 13;
    _EleNum["Si"] = 14;
    _EleNum["P"]  = 15;
    _EleNum["S"]  = 16;
    _EleNum["Cl"] = 17;
    _EleNum["Ar"] = 18;
    _EleNum["K"]  = 19;
    _EleNum["Ca"] = 20;
    _EleNum["Sc"] = 21;
    _EleNum["Ti"] = 22;
    _EleNum["V"]  = 23;
    _EleNum["Cr"] = 24;
    _EleNum["Mn"] = 25;
    _EleNum["Fe"] = 26;
    _EleNum["Co"] = 27;
    _EleNum["Ni"] = 28;
    _EleNum["Cu"] = 29;
    _EleNum["Zn"] = 30;
    _EleNum["Ga"] = 31;
    _EleNum["Ge"] = 32;
    _EleNum["As"] = 33;
    _EleNum["Se"] = 34;
    _EleNum["Br"] = 35;
    _EleNum["Kr"] = 36;
    _EleNum["Rb"] = 37;
    _EleNum["Sr"] = 38;
    _EleNum["Y"]  = 39;
    _EleNum["Zr"] = 40;
    _EleNum["Nb"] = 41;
    _EleNum["Mo"] = 42;
    _EleNum["Tc"] = 43;
    _EleNum["Ru"] = 44;
    _EleNum["Rh"] = 45;
    _EleNum["Pd"] = 46;
    _EleNum["Ag"] = 47;
    _EleNum["Cd"] = 48;
    _EleNum["In"] = 49;
    _EleNum["Sn"] = 50;
    _EleNum["Sb"] = 51;
    _EleNum["Te"] = 52;
    _EleNum["I"]  = 53;
    _EleNum["Xe"] = 54;
    _EleNum["Cs"] = 55;
    _EleNum["Ba"] = 56;
    _EleNum["Hf"] = 72;
    _EleNum["Ta"] = 73;
    _EleNum["W"]  = 74;
    _EleNum["Re"] = 75;
    _EleNum["Os"] = 76;
    _EleNum["Ir"] = 77;
    _EleNum["Pt"] = 78;
    _EleNum["Au"] = 79;
    _EleNum["Hg"] = 80;
    _EleNum["Tl"] = 81;
    _EleNum["Pb"] = 82;
    _EleNum["Bi"] = 83;
    _EleNum["Po"] = 84;
    _EleNum["At"] = 85;
    _EleNum["Rn"] = 86;
  };

  inline void FillEleName() {
    // Nuclear Charges
    _EleName[1]  = "H";
    _EleName[2]  = "He";
    _EleName[3]  = "Li";
    _EleName[4]  = "Be";
    _EleName[5]  = "B";
    _EleName[6]  = "C";
    _EleName[7]  = "N";
    _EleName[8]  = "O";
    _EleName[9]  = "F";
    _EleName[10] = "Ne";
    _EleName[11] = "Na";
    _EleName[12] = "Mg";
    _EleName[13] = "Al";
    _EleName[14] = "Si";
    _EleName[15] = "P";
    _EleName[16] = "S";
    _EleName[17] = "Cl";
    _EleName[18] = "Ar";
    _EleName[19] = "K";
    _EleName[20] = "Ca";
    _EleName[21] = "Sc";
    _EleName[22] = "Ti";
    _EleName[23] = "V";
    _EleName[24] = "Cr";
    _EleName[25] = "Mn";
    _EleName[26] = "Fe";
    _EleName[27] = "Co";
    _EleName[28] = "Ni";
    _EleName[29] = "Cu";
    _EleName[30] = "Zn";
    _EleName[31] = "Ga";
    _EleName[32] = "Ge";
    _EleName[33] = "As";
    _EleName[34] = "Se";
    _EleName[35] = "Br";
    _EleName[36] = "Kr";
    _EleName[37] = "Rb";
    _EleName[38] = "Sr";
    _EleName[39] = "Y";
    _EleName[40] = "Zr";
    _EleName[41] = "Nb";
    _EleName[42] = "Mo";
    _EleName[43] = "Tc";
    _EleName[44] = "Ru";
    _EleName[45] = "Rh";
    _EleName[46] = "Pd";
    _EleName[47] = "Ag";
    _EleName[48] = "Cd";
    _EleName[49] = "In";
    _EleName[50] = "Sn";
    _EleName[51] = "Sb";
    _EleName[52] = "Te";
    _EleName[53] = "I";
    _EleName[54] = "Xe";
    _EleName[55] = "Cs";
    _EleName[56] = "Ba";

    _EleName[72] = "Hf";
    _EleName[73] = "Ta";
    _EleName[74] = "W";
    _EleName[75] = "Re";
    _EleName[76] = "Os";
    _EleName[77] = "Ir";
    _EleName[78] = "Pt";
    _EleName[79] = "Au";
    _EleName[80] = "Hg";
    _EleName[81] = "Tl";
    _EleName[82] = "Pb";
    _EleName[83] = "Bi";
    _EleName[84] = "Po";
    _EleName[85] = "At";
    _EleName[86] = "Rn";
  };

  inline void FillEleShort() {
    // VdW radii in Angstrom as used in MK Gaussian
    _EleShort["HYDROGEN"]   = "H";
    _EleShort["HELIUM"]     = "He";
    _EleShort["LITHIUM"]    = "Li";
    _EleShort["BERYLLIUM"]  = "Be";
    _EleShort["BORON"]      = "B";
    _EleShort["CARBON"]     = "C";
    _EleShort["NITROGEN"]   = "Ni";
    _EleShort["OXYGEN"]     = "O";
    _EleShort["FLUORINE"]   = "Fl";
    _EleShort["NEON"]       = "Ne";
    _EleShort["SODIUM"]     = "Na";
    _EleShort["MAGNESIUM"]  = "Mg";
    _EleShort["ALUMINUM"]   = "Al";
    _EleShort["SILICON"]    = "Si";
    _EleShort["PHOSPHORUS"] = "Ph";
    _EleShort["SULFUR"]     = "S";
    _EleShort["CLORINE"]    = "Cl";
    _EleShort["ARGON"]      = "Ar";
    _EleShort["POTASSIUM"]  = "K";
    _EleShort["CALCIUM"]    = "Ca";
    _EleShort["SCANDDIUM"]  = "Sc";
    _EleShort["TITANIUM"]   = "Ti";
    _EleShort["VANADIUM"]   = "V";
    _EleShort["CHROMIUM"]   = "Cr";
    _EleShort["MANGANESE"]  = "Mn";
    _EleShort["IRON"]       = "Fe";
    _EleShort["COBALT"]     = "Co";
    _EleShort["NICKEL"]     = "Ni";
    _EleShort["COPPER"]     = "Cu";
    _EleShort["ZINC"]       = "Zn";
    _EleShort["GALLIUM"]    = "Ga";
    _EleShort["GERMANIUM"]  = "Ge";
    _EleShort["ARSENIC"]    = "As";
    _EleShort["SELENIUM"]   = "Se";
    _EleShort["BROMINE"]    = "Br";
    _EleShort["KRPTON"]     = "Kr";
    _EleShort["RUBIDIUM"]   = "Rb";
    _EleShort["STRONTIUM"]  = "Sr";
    _EleShort["YTTRIUM"]    = "Y";
    _EleShort["ZIRCONIUM"]  = "Zr";
    _EleShort["NIOBIUM"]    = "Nb";
    _EleShort["MOLYBDENUM"] = "Mo";
    _EleShort["TECHNETIUM"] = "Tc";
    _EleShort["RUTHENIUM"]  = "Ru";
    _EleShort["RHODIUM"]    = "Rh";
    _EleShort["PALLADIUM"]  = "Pd";
    _EleShort["SILVER"]     = "Ag";
    _EleShort["CADMIUM"]    = "Cd";
    _EleShort["INDIUM"]     = "In";
    _EleShort["TIN"]        = "Sn";
    _EleShort["ANTIMONY"]   = "Sb";
    _EleShort["TELLURIUM"]  = "Te";
    _EleShort["IODINE"]     = "I";
    _EleShort["XENON"]      = "Xe";
    _EleShort["CAESIUM"]    = "Cs";
    _EleShort["BARIUM"]     = "Ba";
    _EleShort["HAFNIUM"]    = "Hf";
    _EleShort["TANTALUM"]   = "Ta";
    _EleShort["TUNGSTEN"]   = "W";
    _EleShort["RHENIUM"]    = "Re";
    _EleShort["OSMIUM"]     = "Os";
    _EleShort["IRIDIUM"]    = "Ir";
    _EleShort["PLATINUM"]   = "Pt";
    _EleShort["GOLD"]       = "Au";
    _EleShort["MERCURY"]    = "Hg";
    _EleShort["THALLIUM"]   = "Tl";
    _EleShort["LEAD"]       = "Pd";
    _EleShort["BISMUTH"]    = "Bi";
    _EleShort["PLONIUM"]    = "Po";
    _EleShort["ASTATINE"]   = "At";
    _EleShort["RADON"]      = "Rn";
  };

  inline void FillEleFull() {
    // VdW radii in Angstrom as used in MK Gaussian
    _EleFull["H"]  = "HYDROGEN";
    _EleFull["He"] = "HELIUM";
    _EleFull["Li"] = "LITHIUM";
    _EleFull["Be"] = "BERYLLIUM";
    _EleFull["B"]  = "BORON";
    _EleFull["C"]  = "CARBON";
    _EleFull["N"]  = "NITROGEN";
    _EleFull["O"]  = "OXYGEN";
    _EleFull["F"]  = "FLUORINE";
    _EleFull["Ne"] = "NEON";
    _EleFull["Na"] = "SODIUM";
    _EleFull["Mg"] = "MAGNESIUM";
    _EleFull["Al"] = "ALUMINUM";
    _EleFull["Si"] = "SILICON";
    _EleFull["P"]  = "PHOSPHORUS";
    _EleFull["S"]  = "SULFUR";
    _EleFull["Cl"] = "CHLORINE";
    _EleFull["Ar"] = "ARGON";
    _EleFull["K"]  = "POTASSIUM";
    _EleFull["Ca"] = "CALCIUM";
    _EleFull["Sc"] = "SCANDIUM";
    _EleFull["Ti"] = "TITANIUM";
    _EleFull["V"]  = "VANADIUM";
    _EleFull["Cr"] = "CHROMIUM";
    _EleFull["Mn"] = "MANGANESE";
    _EleFull["Fe"] = "IRON";
    _EleFull["Co"] = "COBALT";
    _EleFull["Ni"] = "NICKEL";
    _EleFull["Cu"] = "COPPER";
    _EleFull["Zn"] = "ZINC";
    _EleFull["Ga"] = "GALLIUM";
    _EleFull["Ge"] = "GERMANIUM";
    _EleFull["As"] = "ARSENIC";
    _EleFull["Se"] = "SELENIUM";
    _EleFull["Br"] = "BROMINE";
    _EleFull["Kr"] = "KRYPTON";
    _EleFull["Rb"] = "RUBIDIUM";
    _EleFull["Sr"] = "STRONTIUM";
    _EleFull["Y"]  = "YTTRIUM";
    _EleFull["Zr"] = "ZIRCONIUM";
    _EleFull["Nb"] = "NIOBIUM";
    _EleFull["Mo"] = "MOLYBDENUM";
    _EleFull["Tc"] = "TECHNETIUM";
    _EleFull["Ru"] = "TUTHENIUM";
    _EleFull["Rh"] = "RHODIUM";
    _EleFull["Pd"] = "PALLADIUM";
    _EleFull["Ag"] = "SILVER";
    _EleFull["Cd"] = "CADMIUM";
    _EleFull["In"] = "INDIUM";
    _EleFull["Sn"] = "TIN";
    _EleFull["Sb"] = "ANTIMONY";
    _EleFull["Te"] = "TELLURIUM";
    _EleFull["I"]  = "IODINE";
    _EleFull["Xe"] = "XENON";
    _EleFull["Cs"] = "CEASIUM";
    _EleFull["Ba"] = "BARIUM";
    _EleFull["Hf"] = "HAFNIUM";
    _EleFull["Ta"] = "TANTALUM";
    _EleFull["W"]  = "TUNGSTEN";
    _EleFull["Re"] = "RHENIUM";
    _EleFull["Os"] = "OSMIUM";
    _EleFull["Ir"] = "IRIDIUM";
    _EleFull["Pt"] = "PLATINUM";
    _EleFull["Au"] = "GOLD";
    _EleFull["Hg"] = "MERCURY";
    _EleFull["Tl"] = "THALLIUM";
    _EleFull["Pb"] = "LEAD";
    _EleFull["Bi"] = "BISMUTH";
    _EleFull["Po"] = "POLONIUM";
    _EleFull["At"] = "ASTATINE";
    _EleFull["Rn"] = "RADON";
  };

  inline void FillVdWMK() {
    // VdW radii in Angstrom as used in MK Gaussian
    _VdWMK["H"]  = 1.2;
    _VdWMK["He"] = 1.2;
    _VdWMK["Li"] = 1.37;
    _VdWMK["Be"] = 1.45;
    _VdWMK["B"]  = 1.45;
    _VdWMK["C"]  = 1.5;
    _VdWMK["N"]  = 1.5;
    _VdWMK["O"]  = 1.4;
    _VdWMK["F"]  = 1.35;
    _VdWMK["Ne"] = 1.3;
    _VdWMK["Na"] = 1.57;
    _VdWMK["Mg"] = 1.36;
    _VdWMK["Al"] = 1.24;
    _VdWMK["Si"] = 1.17;
    _VdWMK["P"]  = 1.8;
    _VdWMK["S"]  = 1.75;
    _VdWMK["Cl"] = 1.7;
    _VdWMK["Ag"] = 2.0;
  };
  
  inline void FillPolarizability(){
    _ElPolarizability["H"] = 0.496e-3;
    _ElPolarizability["C"] = 1.334e-3;
    _ElPolarizability["N"] = 1.073e-3;
    _ElPolarizability["O"] = 0.837e-3;
    _ElPolarizability["S"] = 2.926e-3;
    _ElPolarizability["F"] = 0.440e-3;
    _ElPolarizability["Si"] = 3.962e-3;   // B3LYP/6-311+g(2d,2p)
    _ElPolarizability["Zn"] = 5.962e-3;   // B3LYP/6-311+g(2d,2p)
    _ElPolarizability["Al"] = 5.80e-3;   //[1]P. Fuentealba, “The static dipole polarizability of aluminium atom: discrepancy between theory and experiment,” Chemical physics letters, vol. 397, no. 4, pp. 459–461, 2004.
    };
};
}
}

#endif /* __VOTCA_TOOLS_ELEMENTS_H */
