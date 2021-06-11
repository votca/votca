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

#ifndef VOTCA_TOOLS_ELEMENTS_H
#define VOTCA_TOOLS_ELEMENTS_H

// Standard includes
#include <map>
#include <string>

// Local VOTCA includes
#include "constants.h"

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

  /// Determine if the name is a recognized element symbol or name
  bool isElement(std::string name);

  /// ChelpG is a method for calculating the electrostatic potential outside
  /// of a molecule. The VdWChelpG radius in this case is specific to each
  /// element and defines the radius outside of which the electrostatic
  /// potential is calculated. The electrostatic potential within the radius
  /// is not calculated. For more information see the reference paper CHELPG
  /// paper [Journal of Computational Chemistry 11, 361, 1990]
  double getVdWChelpG(std::string name);

  /// Merz-Singh-Kollman MK method is a similar method for calculating the
  /// electrostatic potential outside of a molecule. Details of the method can
  /// be found here [Singh, U. Chandra (04/01/1984). "An approach to computing
  /// electrostatic charges for molecules". Journal of computational chemistry
  /// 5 (2), 129, 1984]. The VdWMK method will return the radii used to define
  /// where the electrostatic grid starts.
  double getVdWMK(std::string name);

  /// Return the Nuclear charges of each atom. H - 1, He - 2, Na - 3 etc...
  Index getNucCrg(std::string name);

  /// Similar to the Nuclear charges but returns in integer form represents
  /// the id of the atom in the periodic table, the id starts at 1
  Index getEleNum(std::string name);

  /// Returns the mass of each atom in a.u.
  double getMass(std::string name);

  /// Returns the atomic polarisability of atom
  // All polarizabilities in nm**3
  // Isotropic polarizability volume is evaluated from the tensor
  // as (a_xx * a_yy * a_zz )^(1/3) for eigenvalues of the polarizability tensor
  double getPolarizability(std::string name);

  /// Returns the covalent Radii of the atom
  double getCovRad(std::string name, std::string unit);

  /// Provided the element number returns the symbol for the element name
  /// (1) = "H", (2) = "He", ...
  std::string getEleName(Index elenum);

  /// Provided the full element name returns the element symbol
  /// "Hydrogen" = "H", "HELIUM" = "He",...
  std::string getEleShort(std::string elefull);

  /// Is `eleshort` recognized an element symbol i.e. H, C, He, Ne etc
  bool isEleShort(std::string shortname);

  /// Is `elefull` recognized as an element name i.e. Carbon, HYDROGEN, suphur
  bool isEleFull(std::string fullname);

  bool isMassAssociatedWithElement(double mass, double tolerance);

  /// Get the shortened element name given a mass similar in size to one of
  /// the elements. Provided the mass is within the specified tolerance of
  /// the match.
  std::string getEleShortClosestInMass(double mass, double tolerance);

  /// Provided the element symbol returns the element name
  /// "Pb" = "LEAD", "Na" = "SODIUM", ....
  std::string getEleFull(std::string eleshort);

 private:
  // cache variables
  bool filled_VdWChelpG_ = false;
  bool filled_VdWMK_ = false;
  bool filled_NucCrg_ = false;
  bool filled_CovRad_ = false;
  bool filled_Mass_ = false;
  bool filled_EleNum_ = false;
  bool filled_EleName_ = false;
  bool filled_ElPolarizability_ = false;
  bool filled_EleShort_ = false;
  bool filled_EleFull_ = false;

  std::map<std::string, double> VdWChelpG_;
  std::map<std::string, double> VdWMK_;
  std::map<std::string, Index> NucCrg_;
  std::map<std::string, double> CovRad_;
  std::map<std::string, double> Mass_;
  std::map<std::string, Index> EleNum_;
  std::map<Index, std::string> EleName_;

  std::map<std::string, double> ElPolarizability_;

  std::map<std::string, std::string> EleShort_;
  std::map<std::string, std::string> EleFull_;

  /// Finds the element closest in mass and returns the difference as well as
  /// the string of elements short name
  std::pair<std::string, double> findShortNameOfElementClosestInMass(
      double mass);

  void FillMass();
  void FillVdWChelpG();
  void FillNucCrg();
  void FillCovRad();
  void FillEleNum();
  void FillEleName();
  void FillEleShort();
  void FillEleFull();
  void FillVdWMK();
  void FillPolarizability();
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_ELEMENTS_H
