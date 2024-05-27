/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <sstream>
#include <unordered_map>
#include <vector>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/getline.h>

// Local private VOTCA includes
#include "pdbreader.h"

namespace votca {
namespace csg {
using namespace boost;
using namespace std;

bool PDBReader::ReadTopology(string file, Topology &top) {
  topology_ = true;
  top.Cleanup();

  fl_.open(file);
  if (!fl_.is_open()) {
    throw std::ios_base::failure("Error on open topology file: " + file);
  }

  NextFrame(top);
  fl_.close();

  return true;
}

bool PDBReader::Open(const string &file) {

  fl_.open(file);
  if (!fl_.is_open()) {
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  }

  return true;
}

void PDBReader::Close() { fl_.close(); }

bool PDBReader::FirstFrame(Topology &top) {
  topology_ = false;
  return NextFrame(top);
}

bool PDBReader::NextFrame(Topology &top) {

  cout << "The VOTCA pdb reader attempts to follow the official format: "
       << endl;
  cout << "http://www.wwpdb.org/documentation/file-format" << endl;
  cout << endl;
  cout << "NOTE the element symbol is important for associating a mass with "
       << endl;
  cout << "each atom, if no element symbol is specified the atom name will be"
       << endl;
  cout << "used and assumed to be an element symbol." << endl;
  cout << "The pdb file format does not support pseudo atoms such as a CH3 "
       << endl;
  cout << "and thus VOTCA cannot read such atoms in when using a .pdb file."
       << endl;
  cout << "If you need to use pseudo atoms consider using a different file "
       << endl;
  cout << "format such as the lammps data format." << endl;

  string line;
  tools::Elements elements;
  // Two column vector for storing all bonds
  // 1 - id of first atom
  // 2 - id of second atom
  vector<vector<Index>> bond_pairs;
  // Store pointers to every bead
  // WARNING we are assuming in the bead_Eigen::Vector3d that the indices of the
  // beads
  //         correspond to the order in which they are read in. As in the first
  //         bead read in will be at index 0, etc...
  vector<Bead *> bead_vec;
  ////////////////////////////////////////////////////////////////////////////////
  // Read in information from .pdb file
  ////////////////////////////////////////////////////////////////////////////////
  Index bead_count = 0;
  while (tools::getline(fl_, line)) {
    if (tools::wildcmp("CRYST1*", line)) {
      string a, b, c, alpha, beta, gamma;
      try {
        // 1 -  6       Record name    "CRYST1"
        a = string(line, (7 - 1), 9);
        // 7 - 15       Real(9.3)      a           (Angstroms)
        b = string(line, (16 - 1), 9);
        // 16 - 24       Real(9.3)     b           (Angstroms)
        c = string(line, (25 - 1), 9);
        // 25 - 33       Real(9.3)     c           (Angstroms)
        alpha = string(line, (34 - 1), 7);
        // 34 - 40       Real(7.2)     alpha (degrees)
        beta = string(line, (41 - 1), 7);
        // 41 - 47       Real(7.2)     beta        (degrees)
        gamma = string(line, (48 - 1), 7);
        // 48 - 54       Real(7.2)     gamma (degrees)
        // 56 - 66       LString       Space group
        // 67 - 70       Integer       Z value
      } catch (std::out_of_range &) {
        throw std::runtime_error("Misformated pdb file in CRYST1 line");
      }
      boost::algorithm::trim(a);
      boost::algorithm::trim(b);
      boost::algorithm::trim(c);
      boost::algorithm::trim(alpha);
      boost::algorithm::trim(beta);
      boost::algorithm::trim(gamma);
      if ((!tools::wildcmp("90*", alpha)) || (!tools::wildcmp("90*", beta)) ||
          (!tools::wildcmp("90*", gamma))) {
        throw std::runtime_error(
            "Non cubical box in pdb file not implemented, yet!");
      }
      double aa = stod(a) / 10.0;
      double bb = stod(b) / 10.0;
      double cc = stod(c) / 10.0;
      Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
      box.diagonal() = Eigen::Vector3d(aa, bb, cc);
      top.setBox(box);
    }
    // Only read the CONECT keyword if the topology is set too true
    if (topology_ && tools::wildcmp("CONECT*", line)) {
      vector<string> bonded_atms;
      string atm1;
      try {
        // If the CONECT keyword is found then there must be at least
        // two atom identifiers, more than that is optional.

        stringstream ss;
        ss << string(line.substr(6));
        // 1 -  6       Record name    "CONECT"
        // 11 -  7       Real(5)        atm1           (ID)
        ss >> atm1;
        string temp_atm;
        ss >> temp_atm;
        bonded_atms.push_back(temp_atm);
        ss >> temp_atm;
        // Here we have taken a less rigorous approach to the .pdb files
        // we do not care at this point how large the ids of the atoms are
        // they can be greater than 99,999 with this approach.
        while (ss) {
          bonded_atms.push_back(temp_atm);
          ss >> temp_atm;
        }
      } catch (std::out_of_range &) {
        throw std::runtime_error("Misformated pdb file in CONECT line\n" +
                                 line);
      }

      vector<Index> row(2);
      boost::algorithm::trim(atm1);
      Index at1 = boost::lexical_cast<Index>(atm1);
      row.at(0) = at1;

      for (auto &bonded_atm : bonded_atms) {
        Index at2 = boost::lexical_cast<Index>(bonded_atm);
        row.at(1) = at2;
        // Because every bond will be counted twice in a .pdb file
        // we will only add bonds where the id (atm1) is less than the
        // bonded_atm
        if (at1 < at2) {
          bond_pairs.push_back(row);
        }
      }
    }

    if (tools::wildcmp("ATOM*", line) || tools::wildcmp("HETATM*", line)) {

      // according to PDB format
      string x, y, z, resNum, resName, atName;
      string charge, elem_sym;
      // string atNum;
      try {
        /* Some pdb don't include all this, read only what we really need*/
        /* leave this here in case we need more later*/

        // str       ,  "ATOM", "HETATM"
        // string recType    (line,( 1-1),6);
        // Index       , Atom serial number
        // atNum    =    string(line,( 7-1),6);
        // str       , Atom name
        atName = string(line, (13 - 1), 4);
        // char      , Alternate location indicator
        // string atAltLoc   (line,(17-1),1);
        // str       , Residue name
        resName = string(line, (18 - 1), 3);
        // char      , Chain identifier
        // string chainID    (line,(22-1),1);
        // Index       , Residue sequence number
        resNum = string(line, (23 - 1), 4);
        // char      , Code for insertion of res
        // string atICode    (line,(27-1),1);
        // float 8.3 , x
        x = string(line, (31 - 1), 8);
        // float 8.3 , y
        y = string(line, (39 - 1), 8);
        // float 8.3 , z
        z = string(line, (47 - 1), 8);
        // float 6.2 , Occupancy
        // string atOccup    (line,(55-1),6);
        // float 6.2 , Temperature factor
        // string atTFactor  (line,(61-1),6);
        // str       , Segment identifier
        // string segID      (line,(73-1),4);
        // str       , Element symbol
        elem_sym = string(line, (77 - 1), 2);
        // str       , Charge on the atom
        charge = string(line, (79 - 1), 2);
      } catch (std::out_of_range &) {
        string err_msg = "Misformated pdb file in atom line # " +
                         boost::lexical_cast<string>(bead_count) +
                         "\n the correct pdb file format requires 80 "
                         "characters in width (spaces matter). Furthermore, " +
                         "\n to read the topology in from a .pdb file the "
                         "following attributes must be " +
                         "\n specified:                                        "
                         "                        " +
                         "\n Atom Name, Residue Name, Residue Number, x, y, z, "
                         "element symbol"
                         "charge (optional)     \n";
        throw std::runtime_error(err_msg);
      }
      boost::algorithm::trim(atName);
      boost::algorithm::trim(resName);
      boost::algorithm::trim(resNum);
      boost::algorithm::trim(x);
      boost::algorithm::trim(y);
      boost::algorithm::trim(z);
      boost::algorithm::trim(elem_sym);
      boost::algorithm::trim(charge);

      bead_count++;

      Bead *b;
      // Only read the CONECT keyword if the topology is set too true
      if (topology_) {
        Index resnr;
        try {
          resnr = boost::lexical_cast<Index>(resNum);
        } catch (bad_lexical_cast &) {
          throw std::runtime_error(
              "Cannot convert resNum='" + resNum +
              "' to int, that usallly means: misformated pdb file");
        }

        if (resName == "") {
          cout << "WARNING no resname specified, assigning name to: UNK"
               << endl;
          resName = "UNK";
        }

        if (resnr < 1) {
          throw std::runtime_error("Misformated pdb file, resnr has to be > 0");
        }
        // TODO: fix the case that resnr is not in ascending order
        if (resnr > top.ResidueCount()) {
          while ((resnr - 1) > top.ResidueCount()) {  // pdb resnr should start
                                                      // with 1 but accept
                                                      // sloppy files

            // create dummy residue, hopefully it will never show

            top.CreateResidue(resName);
            cout << "Warning: residue numbers not continuous, create dummy "
                    "residue with residue number "
                 << top.ResidueCount() << endl;
          }
          top.CreateResidue(resName);
        }
        // This is not correct, but still better than no type at all!
        if (!top.BeadTypeExist(atName)) {
          top.RegisterBeadType(atName);
        }

        // Determine if the charge has been provided in the .pdb file or if we
        // will be assuming it is 0
        double ch = 0;
        if (charge != "") {
          ch = stod(charge);
        } else {
          cout << "WARNING no charge was specified for " << endl;
          cout << line << endl;
          cout << "Assuming a charge of 0" << endl;
        }

        if (elem_sym == "") {
          cout << "WARNING no element was specified, assuming atom name is "
               << endl;
          cout << "an element symbol: " << atName << endl;
          if (elements.isElement(atName)) {
            elem_sym = atName;
          } else {
            throw std::runtime_error(
                "Atom name is not an element symbol, so substitution fails, "
                "the element is needed in order to resolve the mass.");
          }
        }

        // CreateBead takes 6 parameters in the following order
        // 1 - symmetry of the bead
        // 2 - name of the bead     (string)
        // 3 - bead type            (BeadType *)
        // 4 - residue number       (Index)
        // 5 - mass                 (double)
        // 6 - charge               (double)
        //
        // res -1 as internal number starts with 0
        b = top.CreateBead(Bead::spherical, atName, atName, resnr - 1,
                           elements.getMass(elem_sym), ch);
      } else {
        b = top.getBead(bead_count - 1);
      }
      // convert to nm from A
      b->setPos(
          Eigen::Vector3d(stod(x) / 10.0, stod(y) / 10.0, stod(z) / 10.0));

      bead_vec.push_back(b);
    }

    if ((line == "ENDMDL") || (line == "END") || (fl_.eof())) {
      break;
    }
  }

  if (!topology_ && (bead_count > 0) && bead_count != top.BeadCount()) {
    throw std::runtime_error(
        "number of beads in topology and trajectory differ");
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Sort data and determine atom structure, connect with top (molecules, bonds)
  ////////////////////////////////////////////////////////////////////////////////

  // Extra processing is done if the file is a topology file, in which case the
  // atoms must be sorted into molecules and the bond interactions recorded
  if (topology_) {
    // Now we need to add the bond pairs
    // WARNING We are assuming the atom ids are contiguous with no gaps

    // First Index  - is the index of the atom
    // Second Index - is the index of the molecule
    map<Index, Index> atm_molecule;

    // First Index  - is the index of the molecule
    // list<Index>  - is a list of the atoms in the molecule
    unordered_map<Index, list<Index>> molecule_atms;

    // Keep track of the number of molecules we have created through an index
    Index mol_index = 0;

    // Cycle through all bonds
    for (auto &bond_pair : bond_pairs) {

      Index atm_id1 = bond_pair.at(0);
      Index atm_id2 = bond_pair.at(1);
      // Check to see if either atm referred to in the bond is already
      // attached to a molecule
      auto mol_iter1 = atm_molecule.find(atm_id1);
      auto mol_iter2 = atm_molecule.find(atm_id2);

      // This means neither atom is attached to a molecule
      if (mol_iter1 == atm_molecule.end() && mol_iter2 == atm_molecule.end()) {
        // We are going to create a new row for a new molecule
        list<Index> atms_in_mol;
        atms_in_mol.push_back(atm_id1);
        atms_in_mol.push_back(atm_id2);
        molecule_atms[mol_index] = atms_in_mol;
        // Associate atm1 and atm2 with the molecule index
        atm_molecule[atm_id1] = mol_index;
        atm_molecule[atm_id2] = mol_index;
        // Increment the molecule index
        mol_index++;

        // This means only atm2 is attached to a molecule
      } else if (mol_iter1 == atm_molecule.end()) {
        // Add atm1 to the molecule that contains atm2
        molecule_atms[mol_iter2->second].push_back(atm_id1);
        // Associate atm1 with the molecule it is now part of
        atm_molecule[atm_id1] = mol_iter2->second;

        // This means only atm1 is attached to a molecule
      } else if (mol_iter2 == atm_molecule.end()) {
        // Add atm2 to the molecule that contains atm1
        molecule_atms[mol_iter1->second].push_back(atm_id2);
        // Associate atm1 with the molecule it is now part of
        atm_molecule[atm_id2] = mol_iter1->second;

      } else if (mol_iter1 != mol_iter2) {
        // This means both atm1 and atm2 are attached to a molecule
        // But if they are already attached to the same molecule there is
        // nothing else to be done.
        Index chosen_mol;
        Index obsolete_mol;
        // We will merge the atms to the molecule with the smallest index
        if (mol_iter1->second < mol_iter2->second) {
          chosen_mol = mol_iter1->second;
          obsolete_mol = mol_iter2->second;
        } else {
          chosen_mol = mol_iter2->second;
          obsolete_mol = mol_iter1->second;
        }

        // Now we will proceed to cycle through the atms that were in the now
        // obsolete molecule and make sure they are pointing to the new molecule
        for (Index &atm_temp : molecule_atms[obsolete_mol]) {

          atm_molecule[atm_temp] = chosen_mol;
        }

        // Splicing will remove atoms from the now obsolete molecule and place
        // them on the chosen molecule.
        molecule_atms[chosen_mol].splice(molecule_atms[chosen_mol].end(),
                                         molecule_atms[obsolete_mol]);

        // Finally we will clear out the record of the obsolete molecule
        molecule_atms.erase(obsolete_mol);
      }
    }
#ifndef NDEBUG
    cerr << "Consistency check for pdbreader" << endl;
    Index i = 0;
    for (auto lis = molecule_atms.begin(); lis != molecule_atms.end(); lis++) {
      cerr << "Molecule " << i << endl;
      cerr << "Atoms: ";
      for (auto atm_ind = lis->second.begin(); atm_ind != lis->second.end();
           atm_ind++) {
        cerr << *atm_ind << " ";
      }
      cerr << endl;
      i++;
    }
    cerr << endl;
#endif
    // Now that we know which interactions belong to which molecules we can:
    // 1 Add the molecules
    // 2 Add the bond interactions

    // Molecule map
    // First Index - is the index of the molecule
    // Molecule* - is a pointer to the Molecule object
    map<Index, Molecule *> mol_map;

    // Used to reindex the molecules so that they start at 0 and progress
    // with out gaps in their ids.
    // First Index  - is the index of the old molecule
    // Second Index - is the new index
    map<Index, Index> mol_reInd_map;

    Index ind = 0;
    for (auto mol = molecule_atms.begin(); mol != molecule_atms.end(); mol++) {

      string mol_name = "PDB Molecule " + boost::lexical_cast<string>(ind);

      Molecule *mi = top.CreateMolecule(mol_name);
      mol_map[mol->first] = mi;
      mol_reInd_map[mol->first] = ind;

      // Add all the atoms to the appropriate molecule object
      list<Index> atm_list = molecule_atms[mol->first];
      for (Index &atm_temp : atm_list) {

        string residuename = "DUM";
        mi->AddBead(bead_vec.at(atm_temp - 1), residuename);
      }
      ind++;
    }

    Index bond_indx = 0;
    // Cyle through the bonds and add them to the appropriate molecule
    for (auto &bond_pair : bond_pairs) {

      Index atm_id1 = bond_pair.at(0);
      Index atm_id2 = bond_pair.at(1);
      // Should be able to just look at one of the atoms the bond is attached
      // too because the other will also be attached to the same molecule.
      Index mol_ind = atm_molecule[atm_id1];
      Molecule *mi = mol_map[mol_ind];
      // Grab the id of the bead associated with the atom
      // It may be the case that the atom id's and bead id's are different
      Index bead_id1 = bead_vec.at(atm_id1 - 1)->getId();
      Index bead_id2 = bead_vec.at(atm_id2 - 1)->getId();
      Interaction *ic = new IBond(bead_id1, bead_id2);
      ic->setGroup("BONDS");
      ic->setIndex(bond_indx);
      bond_indx++;
      ic->setMolecule(mol_reInd_map[mol_ind]);
      top.AddBondedInteraction(ic);
      mi->AddInteraction(ic);
    }

    // Finally we want to build an exclusion matrix
    top.RebuildExclusions();
  }

  return !fl_.eof();
}
}  // namespace csg
}  // namespace votca
