/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
using namespace votca::tools;

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
  // .at(0) - element symbol
  // .at(1) - atom name may be the same as the element symbol
  std::map<int, std::vector<std::string>> atomtypes_;

  // String is the type .e.g. "atom","bond" etc
  // int is the number of different types
  std::map<string, int> numberOfDifferentTypes_;

  // String is the type .e.g. "atom", "bond" etc
  // int is the number of them
  std::map<string, int> numberOf_;

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

  void InitializeAtomTypes_();
  void ReadBox_(std::vector<std::string> fields, Topology &top);
  void SortIntoDataGroup_(std::string tag);
  void ReadNumTypes_(std::vector<string> fields, string type);

  std::vector<std::string> TrimCommentsFrom_(std::vector<std::string> fields);

  void ReadNumOfAtoms_(std::vector<std::string> fields, Topology &top);
  void ReadNumOfBonds_(std::vector<std::string> fields);
  void ReadNumOfAngles_(std::vector<std::string> fields);
  void ReadNumOfDihedrals_(std::vector<std::string> fields);
  void ReadNumOfImpropers_(std::vector<std::string> fields);

  void ReadAtoms_(Topology &top);
  void ReadBonds_(Topology &top);
  void ReadAngles_(Topology &top);
  void ReadDihedrals_(Topology &top);
  void ReadImpropers_(Topology &top);
};
}
}

#endif // _VOTCA_CSG_LAMMPSDATAREADER_H
