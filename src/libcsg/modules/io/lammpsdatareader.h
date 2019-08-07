/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_CSG_LAMMPSDATAREADER_H
#define _VOTCA_CSG_LAMMPSDATAREADER_H

#include <fstream>
#include <iostream>
#include <string>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>

namespace votca {
namespace csg {

class Molecule;

/**
    \brief class for reading lammps data files

    This class provides the TrajectoryReader + Topology reader interface
    for lammps data files

*/
class LAMMPSDataReader : public TrajectoryReader, public TopologyReader {
 public:
  LAMMPSDataReader() {}
  ~LAMMPSDataReader() {}

  /// open, read and close topology file
  bool ReadTopology(std::string file, Topology &top);

  /// open a trajectory file
  bool Open(const std::string &file);
  /// read in the first frame of trajectory file
  bool FirstFrame(Topology &top);
  /// read in the next frame of trajectory file
  bool NextFrame(Topology &top);
  /// close the topology file
  void Close();

 private:
  std::ifstream fl_;
  std::string fname_;
  bool topology_;

  std::map<std::string, std::vector<std::vector<std::string>>> data_;

  // int - atom type starting index of 0
  // .at(0) - element symbol or bead
  // .at(1) - atom name may be the same as the element symbol or bead depending
  //          on if there is more than one atom type for a given element
  std::map<int, std::array<std::string, 2>> atomtypes_;

  // String is the type .e.g. "atom","bond" etc
  // int is the number of different types
  std::map<std::string, int> numberOfDifferentTypes_;

  // String is the type .e.g. "atom", "bond" etc
  // int is the number of them
  std::map<std::string, int> numberOf_;

  // First int is the molecule id
  // Second int is the molecule ptr
  std::map<int, Molecule *> molecules_;

  // First int is the atom id second the molecule id
  std::map<int, int> atomIdToMoleculeId_;

  // First int is the atom id second the atom index
  std::map<int, int> atomIdToIndex_;

  bool MatchOneFieldLabel_(std::vector<std::string> fields, Topology &top);
  bool MatchTwoFieldLabels_(std::vector<std::string> fields, Topology &top);
  bool MatchThreeFieldLabels_(std::vector<std::string> fields, Topology &top);
  bool MatchFourFieldLabels_(std::vector<std::string> fields, Topology &top);
  bool MatchFieldsTimeStepLabel_(std::vector<std::string> fields,
                                 Topology &top);

  void ReadBox_(std::vector<std::string> fields, Topology &top);
  void SortIntoDataGroup_(std::string tag);
  void ReadNumTypes_(std::vector<std::string> fields, std::string type);

  void ReadNumOfAtoms_(std::vector<std::string> fields, Topology &top);
  void ReadNumOfBonds_(std::vector<std::string> fields);
  void ReadNumOfAngles_(std::vector<std::string> fields);
  void ReadNumOfDihedrals_(std::vector<std::string> fields);
  void ReadNumOfImpropers_(std::vector<std::string> fields);

  void ReadAtoms_(Topology &top);
  void ReadBonds_(Topology &top);
  void ReadAngles_(Topology &top);
  void ReadDihedrals_(Topology &top);
  void SkipImpropers_();

  void RenameMolecule(Molecule &mol) const;

  enum lammps_format {
    style_angle_bond_molecule = 0,
    style_atomic = 1,
    style_full = 2
  };
  lammps_format determineDataFileFormat_(std::string line);

  /**
   * \brief Determines atom and bead types based on masses in lammps files
   *
   * The purpose of this function is to take lammps output where there are more
   * than a single atom type of the same element. For instance there may be 4
   * atom types with mass of 12.01. Well this means that they are all carbon but
   * are treated differently in lammps. It makes sense to keep track of this. If
   * a mass cannot be associated with an element we will assume it is pseudo
   * atom or course grained watom which we will represent as a bead. So
   * when creating the atom names we will take this into account. So say we have
   *the following masses in the lammps .data file:
   *
   * Masses
   *
   * 1 1.0
   * 2 12.01
   * 3 12.01
   * 4 16.0
   * 5 12.01
   * 6 15.2
   * 7 12.8
   * 8 15.2
   *
   * Then we would translate this to the following atom names
   * 1 H
   * 2 C1
   * 3 C2
   * 4 O
   * 5 C3
   * 6 Bead1 BeadType 1
   * 7 Bead2 BeadType
   * 8 Bead1 BeadType 2
   *
   * Note that we do not append a number if it is singular, in such cases the
   * element and the atom name is the same.
   **/
  void InitializeAtomAndBeadTypes_();
  std::map<std::string, double> determineBaseNameAssociatedWithMass_();
  std::map<std::string, int> determineAtomAndBeadCountBasedOnMass_(
      std::map<std::string, double> baseNamesAndMasses);
};
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_LAMMPSDATAREADER_H
