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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE lammpdumpstrajectoryreaderwriter_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/elements.h>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

// Check if file exists
bool fexists(const string filename) {
  std::ifstream ifile(filename);
  return (bool)ifile;
}

BOOST_AUTO_TEST_SUITE(lammpsdumpreaderwriter_test)

/**
* \brief Test the trajectory reader
*
* This test is designed to test the trajectory reader this is done by
* creating a small lammps dump file. A topology object is created with
* some default values. The file is then read in with the
 * trajectory reader and the values in the top object are then examined
* to ensure they no longer represent the default state but the values
* from the file.
*/
BOOST_AUTO_TEST_CASE(test_trajectoryreader) {

  // Create a .dump file with (2-bonded thiophene monomers)
  // and read from it. Create a topology object with the same
  // molecule to enable the ability to read in the trajectory
  // file
  string lammpsdumpfilename = "test_thiophene.dump";
  if (fexists(lammpsdumpfilename)) {
    remove(lammpsdumpfilename.c_str());
  }
  ofstream outfile(lammpsdumpfilename);

  // Atom Coordinates
  vector<vector<double>> atom_xyz_file = {
      {3.57166300, -2.91232800, 0.62799100},
      {2.70277100, -1.74647700, 0.72415900},
      {1.44813500, -2.02583600, 0.30856000},
      {2.94693900, -3.99290500, 0.16108100},
      {4.61160400, -2.88232000, 0.91932700},
      {1.28482000, -3.69221600, -0.19709400},
      {3.37487400, -4.97110400, 0.00494100},
      {3.08267998, -0.80674791, 1.09707295},
      {0.24459700, -1.15927300, 0.24162000},
      {-0.59908296, -0.82091411, -1.24794356},
      {-0.36240496, -0.56711816, 1.28476993},
      {-0.01775039, -0.64275306, 2.30554425},
      {-1.54246635, 0.20136386, 0.91593361},
      {-1.76480547, 0.13871054, -0.40089364},
      {-2.56388268, 0.61198783, -0.94726998},
      {-2.12173485, 0.73285944, 1.65650560}};

  // Forces
  vector<vector<double>> atom_forces_file = {
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}};

  outfile << "ITEM: TIMESTEP" << endl;
  outfile << "1" << endl;
  outfile << "ITEM: NUMBER OF ATOMS" << endl;
  outfile << "16" << endl;
  outfile << "ITEM: BOX BOUNDS pp pp pp" << endl;
  outfile << "0 0.000000" << endl;
  outfile << "0 10.000000" << endl;
  outfile << "0 10.000000" << endl;
  outfile << "ITEM: ATOMS id type x y z fx fy fz" << endl;
  outfile << "1 0 3.571663 -2.912328 0.627991 0.100000 0.100000 0.100000"
          << endl;
  outfile << "2 0 2.702771 -1.746477 0.724159 0.100000 0.100000 0.100000"
          << endl;
  outfile << "3 0 1.448135 -2.025836 0.308560 0.100000 0.100000 0.100000"
          << endl;
  outfile << "4 0 2.946939 -3.992905 0.161081 0.100000 0.100000 0.100000"
          << endl;
  outfile << "5 1 4.611604 -2.882320 0.919327 0.100000 0.100000 0.100000"
          << endl;
  outfile << "6 2 1.284820 -3.692216 -0.197094 0.100000 0.100000 0.100000"
          << endl;
  outfile << "7 1 3.374874 -4.971104 0.004941 0.100000 0.100000 0.100000"
          << endl;
  outfile << "8 1 3.082680 -0.806748 1.097073 0.100000 0.100000 0.100000"
          << endl;
  outfile << "9 0 0.244597 -1.159273 0.241620 0.100000 0.100000 0.100000"
          << endl;
  outfile << "10 2 -0.599083 -0.820914 -1.247944 0.100000 0.100000 0.100000"
          << endl;
  outfile << "11 0 -0.362405 -0.567118 1.284770 0.100000 0.100000 0.100000"
          << endl;
  outfile << "12 1 -0.017750 -0.642753 2.305544 0.100000 0.100000 0.100000"
          << endl;
  outfile << "13 0 -1.542466 0.201364 0.915934 0.100000 0.100000 0.100000"
          << endl;
  outfile << "14 0 -1.764805 0.138711 -0.400894 0.100000 0.100000 0.100000"
          << endl;
  outfile << "15 1 -2.563883 0.611988 -0.947270 0.100000 0.100000 0.100000"
          << endl;
  outfile << "16 1 -2.121735 0.732859 1.656506 0.100000 0.100000 0.100000"
          << endl;

  outfile.close();

  Topology top;

  // Make square box
  matrix box;
  // 10  0  0
  //  0 10  0
  //  0  0 10
  box.ZeroMatrix();
  box.set(0, 0, 10);
  box.set(1, 1, 10);
  box.set(2, 2, 10);

  top.setBox(box);
  auto boxType = top.getBoxType();
  BOOST_CHECK_EQUAL(boxType, BoundaryCondition::typeOrthorhombic);
  OrthorhombicBox ortho_box;
  ortho_box.setBox(box);
  top.setStep(1);

  // atom type - this may or may not be the same as the element name, e.g.
  // different atom types depending on the forcefield and bond structure.
  vector<string> atom_types = {"C", "C", "C", "C", "H", "S", "H", "H",
                               "C", "S", "C", "H", "C", "C", "H", "H"};

  // Atom Coordinates
  vector<vector<double>> atom_xyz = {{3.57166300, -1.91232800, 1.62799100},
                                     {2.70277100, -0.74647700, 1.72415900},
                                     {1.44813500, -1.02583600, 1.30856000},
                                     {2.94693900, -2.99290500, 1.16108100},
                                     {4.61160400, -1.88232000, 1.91932700},
                                     {1.28482000, -2.69221600, 0.19709400},
                                     {3.37487400, -3.97110400, 1.00494100},
                                     {3.08267998, 0.80674791, 2.09707295},
                                     {0.24459700, -0.15927300, 1.24162000},
                                     {0.59908296, 0.82091411, -0.24794356},
                                     {0.36240496, 0.56711816, 2.28476993},
                                     {0.01775039, 0.64275306, 3.30554425},
                                     {-0.54246635, 1.20136386, 1.91593361},
                                     {-0.76480547, 1.13871054, 0.40089364},
                                     {-1.56388268, 1.61198783, 0.94726998},
                                     {-1.12173485, 1.73285944, 2.65650560}};

  // Atom Velocities
  vector<vector<double>> atom_vel = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  // Forces
  vector<vector<double>> atom_forces = {
      {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2},
      {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2},
      {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2},
      {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.2, 0.2}};

  Elements elements;
  int residue_num = 1;
  double charge = 0.0;
  byte_t symmetry = 1;

  for (size_t ind = 0; ind < atom_types.size(); ++ind) {
    BeadType * type = top.GetOrCreateBeadType(atom_types.at(ind));
    Bead *b = top.CreateBead(symmetry, atom_types.at(ind), type, residue_num,
                             elements.getMass(atom_types.at(ind)), charge);

    vec xyz(atom_xyz.at(ind).at(0), atom_xyz.at(ind).at(1),
            atom_xyz.at(ind).at(2));
    b->setPos(xyz);

    vec xyz_vel(atom_vel.at(ind).at(0), atom_vel.at(ind).at(1),
                atom_vel.at(ind).at(2));
    b->setVel(xyz_vel);

    vec xyz_forces(atom_forces.at(ind).at(0), atom_forces.at(ind).at(1),
                   atom_forces.at(ind).at(2));
    b->setF(xyz_forces);
  }

  TrajectoryReader::RegisterPlugins();
  TrajectoryReader *reader;
  reader = TrjReaderFactory().Create(lammpsdumpfilename);
  reader->Open(lammpsdumpfilename);
  reader->FirstFrame(top);
  reader->Close();

  for (size_t ind = 0; ind < atom_types.size(); ++ind) {
    Bead *b = top.getBead(ind);
    BOOST_CHECK_CLOSE(b->Pos().x(), atom_xyz_file.at(ind).at(0), 0.01);
    BOOST_CHECK_CLOSE(b->Pos().y(), atom_xyz_file.at(ind).at(1), 0.01);
    BOOST_CHECK_CLOSE(b->Pos().z(), atom_xyz_file.at(ind).at(2), 0.01);
    BOOST_CHECK_CLOSE(b->F().x(), atom_forces_file.at(ind).at(0), 0.01);
    BOOST_CHECK_CLOSE(b->F().y(), atom_forces_file.at(ind).at(1), 0.01);
    BOOST_CHECK_CLOSE(b->F().z(), atom_forces_file.at(ind).at(2), 0.01);
  }
}

/**
 * \brief Testing trajectory writer
 *
 * This test first creates a topology object and assigns default values to
 * it. It then writes the topology info to a lammps dump file. The dump
       * file is then read into the topology file and the values are compared.
 */
BOOST_AUTO_TEST_CASE(test_trajectorywriter) {

  // Create a topology object with a simple system (2-bonded thiophene monomers)
  // and write it to a lammps dump file
  Topology top;

  // Make square box
  matrix box;
  // 10  0  0
  //  0 10  0
  //  0  0 10
  box.ZeroMatrix();
  box.set(0, 0, 10);
  box.set(1, 1, 10);
  box.set(2, 2, 10);

  top.setBox(box);
  auto boxType = top.getBoxType();
  BOOST_CHECK_EQUAL(boxType, BoundaryCondition::typeOrthorhombic);

  OrthorhombicBox ortho_box;
  ortho_box.setBox(box);
  top.setStep(1);

  // atom type - this may or may not be the same as the element name, e.g.
  // different atom types depending on the forcefield and bond structure.
  vector<string> atom_types = {"C", "C", "C", "C", "H", "S", "H", "H",
                               "C", "S", "C", "H", "C", "C", "H", "H"};

  // Atom Coordinates
  vector<vector<double>> atom_xyz = {{3.57166300, -2.91232800, 0.62799100},
                                     {2.70277100, -1.74647700, 0.72415900},
                                     {1.44813500, -2.02583600, 0.30856000},
                                     {2.94693900, -3.99290500, 0.16108100},
                                     {4.61160400, -2.88232000, 0.91932700},
                                     {1.28482000, -3.69221600, -0.19709400},
                                     {3.37487400, -4.97110400, 0.00494100},
                                     {3.08267998, -0.80674791, 1.09707295},
                                     {0.24459700, -1.15927300, 0.24162000},
                                     {-0.59908296, -0.82091411, -1.24794356},
                                     {-0.36240496, -0.56711816, 1.28476993},
                                     {-0.01775039, -0.64275306, 2.30554425},
                                     {-1.54246635, 0.20136386, 0.91593361},
                                     {-1.76480547, 0.13871054, -0.40089364},
                                     {-2.56388268, 0.61198783, -0.94726998},
                                     {-2.12173485, 0.73285944, 1.65650560}};

  // Atom Velocities
  vector<vector<double>> atom_vel = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  // Forces
  vector<vector<double>> atom_forces = {
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}};

  Elements elements;
  int residue_num = 1;
  double charge = 0.0;
  byte_t symmetry = 1;

  for (size_t ind = 0; ind < atom_types.size(); ++ind) {
    BeadType * type = top.GetOrCreateBeadType(atom_types.at(ind));
    Bead *b = top.CreateBead(symmetry, atom_types.at(ind), type, residue_num,
                             elements.getMass(atom_types.at(ind)), charge);

    vec xyz(atom_xyz.at(ind).at(0), atom_xyz.at(ind).at(1),
            atom_xyz.at(ind).at(2));
    b->setPos(xyz);

    vec xyz_vel(atom_vel.at(ind).at(0), atom_vel.at(ind).at(1),
                atom_vel.at(ind).at(2));
    b->setVel(xyz_vel);

    vec xyz_forces(atom_forces.at(ind).at(0), atom_forces.at(ind).at(1),
                   atom_forces.at(ind).at(2));
    b->setF(xyz_forces);
  }
  top.SetHasForce(true);
  top.SetHasVel(true);

  string lammpsDumpFileName = "test_thiophene.dump";
  if (fexists(lammpsDumpFileName)) {
    remove(lammpsDumpFileName.c_str());
  }

  // Write the topology object to a file
  TrajectoryWriter::RegisterPlugins();
  TrajectoryWriter *writer;
  writer = TrjWriterFactory().Create(lammpsDumpFileName);
  writer->Open(lammpsDumpFileName);
  writer->Write(&top);
  writer->Close();

  // Read the .dump file and ensure the information is correct
  TrajectoryReader::RegisterPlugins();
  TrajectoryReader *reader;
  reader = TrjReaderFactory().Create(lammpsDumpFileName);
  reader->Open(lammpsDumpFileName);
  reader->FirstFrame(top);
  reader->Close();

  for (size_t ind = 0; ind < atom_types.size(); ++ind) {
    Bead *b = top.getBead(ind);
    BOOST_CHECK_CLOSE(b->Pos().x(), atom_xyz.at(ind).at(0), 0.01);
    BOOST_CHECK_CLOSE(b->Pos().y(), atom_xyz.at(ind).at(1), 0.01);
    BOOST_CHECK_CLOSE(b->Pos().z(), atom_xyz.at(ind).at(2), 0.01);
    BOOST_CHECK_CLOSE(b->F().x(), atom_forces.at(ind).at(0), 0.01);
    BOOST_CHECK_CLOSE(b->F().y(), atom_forces.at(ind).at(1), 0.01);
    BOOST_CHECK_CLOSE(b->F().z(), atom_forces.at(ind).at(2), 0.01);
  }
}

BOOST_AUTO_TEST_SUITE_END()
