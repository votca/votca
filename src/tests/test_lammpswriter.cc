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

#define BOOST_TEST_MODULE pdbreader_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>
#include <string>
#include <votca/tools/matrix.h>
#include <votca/csg/bead.h>
#include <votca/csg/topologywriter.h>
#include <votca/tools/elements.h>


using namespace std;
using namespace votca::csg;
using namepsace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
	v *= pow(10, p);
	v = round(v);
	v /= pow(10, p);
	return v;
}

BOOST_AUTO_TEST_SUITE(pdbreader_test)

	BOOST_AUTO_TEST_CASE(test_topologywriter) {

		// Create a topology object with a simple system (2-bonded thiophene monomers)
		// and write it to a lammps dump file

   	// Make square box
		matrix box;
		// 10  0  0
		//  0 10  0
		//  0  0 10 
		box.ZeroMatrix();
		box.set(1,1,10);
		box.set(2,2,10);
		box.set(3,3,10);

		top.setBox(box);
		auto boxType = top.getBoxType();

		BoundaryCondition bc(box);
		top.setStep(1);

		// atom type - this may or may not be the same as the element name, e.g. 
		// different atom types depending on the forcefield and bond structure.
		vector<string> atom_types = {
			"C","C","C","C","H","S","H","H",
			"C","S","C","H","C","C","H","H"};

		// Atom Coordinates
		vector<vector<double>> atom_xyz = {
			{3.57166300,   -2.91232800,    0.62799100},
			{2.70277100,   -1.74647700,    0.72415900},
			{1.44813500,   -2.02583600,    0.30856000},
			{2.94693900,   -3.99290500,    0.16108100},
			{4.61160400,   -2.88232000,    0.91932700},
			{1.28482000,   -3.69221600,   -0.19709400},
			{3.37487400,   -4.97110400,    0.00494100},
			{3.08267998,   -0.80674791,    1.09707295},
			{0.24459700,   -1.15927300,    0.24162000},
			{-0.59908296,   -0.82091411,   -1.24794356},
			{-0.36240496,   -0.56711816,    1.28476993},
			{-0.01775039,   -0.64275306,    2.30554425},
			{-1.54246635,    0.20136386,    0.91593361},
			{-1.76480547,    0.13871054,   -0.40089364},
			{-2.56388268,    0.61198783,   -0.94726998},
			{-2.12173485,    0.73285944,    1.65650560}};

		// Atom Velocities
		vector<vector<double>> atom_vel = {
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0 }};

		// Forces
		vector<vector<double>> atom_forces = {
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 },
			{ 0.1, 0.1, 0.1 }};

		for(auto ind=0; ind<atom_types.size();++ind){
			// Atom id starts at 1 but the internals start at 0
			Bead *b = top.getBead(ind);
			b->Pos().x() = atom_xyz.at(ind).at(0);
			b->Pos().y() = atom_xyz.at(ind).at(1);
			b->Pos().z() = atom_xyz.at(ind).at(2);
			b->Vel().x() = atom_vel.at(ind).at(0);
			b->Vel().y() = atom_vel.at(ind).at(1);
			b->Vel().z() = atom_vel.at(ind).at(2);
			b->F().x() = atom_forces.at(ind).at(0);
			b->F().y() = atom_forces.at(ind).at(1);
			b->F().z() = atom_forces.at(ind).at(2);
			BeadType *type = top.GetOrCreateBeadType(atom_types.at(ind));
			b->setType(type);
		}

		TopologyWriter::RegisterPlugins();
		TopologyWriter* writer;
		string lammpsDumpFileName = "test_thiophene.dump";
		writer = TopWriterFactory().Create(lammpsDumpFileName);
		writer->WriteTopology(lammpsDumpFileName, top);
	}

BOOST_AUTO_TEST_SUITE_END()
