/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp> 
#include <votca/tools/getline.h>

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "dlpolytopologyreader.h"

namespace votca { namespace csg {

using namespace std;


string DLPOLYTopologyReader::_NextKeyline(ifstream &fs, const char* wspace)
{
  string line;
  size_t i_nws=0;

  do {
    getline(fs, line);

    if (fs.eof())
      throw std::runtime_error("Error: unexpected end of dlpoly file '" + _fname + "'");

    i_nws = line.find_first_not_of(wspace);
  } while ( line.substr(i_nws,1) == "#" || line.substr(i_nws,1) == ";" );

  return line.substr(i_nws,line.size()-i_nws);
}

string DLPOLYTopologyReader::_NextKeyInt(ifstream &fs, const char* wspace, const string &word, int &ival)
{
  stringstream sl(_NextKeyline(fs,wspace));
  string line,sval;

  sl >> line; // allow user not to bother about the case
  boost::to_upper(line);

  if( line.substr(0,word.size()) != word )
    throw std::runtime_error("Error: unexpected line in dlpoly file " + _fname + ", expected '" + word + "' but got '" + line + "'");

  sl >> sval;

  size_t i_num = sval.find_first_of("0123456789"); // assume integer number straight after the only keyword

  if( i_num>0 ) 
    throw std::runtime_error("Error: missing integer number in directive '" + line + "' in topology file '"+ _fname +"'");

  ival = boost::lexical_cast<int>(sval);

  return sl.str();
}

bool DLPOLYTopologyReader::_isKeyInt(const string &line, const char* wspace, const string &word, int &ival)
{
  Tokenizer tok(line,wspace); // split directives consisting of a few words - the keyword is the first one!
  vector<string> fields;
  tok.ToVector(fields);

  ival = 0;

  if( fields.size() < 2 ) return false;
  //throw std::runtime_error("Error: missing integer number in directive '" + line + "' in topology file '"+ _fname +"'");

  boost::to_upper(fields[0]);

  if (fields[0].substr(0,word.size()) != word )
    throw std::runtime_error("Error: unexpected directive in dlpoly file '" + _fname + "', expected keyword '"+ word +"' but got '" + fields[0] + "'");

  size_t i_num=string::npos;

  int i=1;
  do { // find integer number in the field with the lowest index (closest to the keyword)
    i_num = fields[i++].find_first_of("0123456789");
  } while ( i_num>0 && i<fields.size() );

  if( i_num>0 ) return false;
  //throw std::runtime_error("Error: missing integer number in directive '" + line + "' in topology file '"+ _fname +"'");

  ival = boost::lexical_cast<int>(fields[i-1]);

  return true;
}


bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top)
{
  //#define WhiteSpace " \t"

  const char *WhiteSpace=" \t";

    int  mavecs=0;
    int  mpbct=0;
    int  matoms=0;
    int  natoms=0;
    bool hasVs=false;
    bool hasFs=false;

    BoundaryCondition::eBoxtype pbc_type=BoundaryCondition::typeAuto;

    std::ifstream fl;
    boost::filesystem::path filepath(file.c_str());

    string line;

    // TODO: fix residue naming / assignment - DL_POLY has no means to recognise residues!
    Residue *res = top.CreateResidue("no");

    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        _fname="FIELD"; // DL_POLY uses fixed file names in current/working directory
      } else {
	_fname=filepath.parent_path().string() + "/FIELD";
      }
    } else {
      _fname=file;
    }

    fl.open(_fname.c_str());

    if (!fl.is_open()){
      throw std::runtime_error("Error on opening dlpoly file '" + _fname + "'");
    } else {

      line = _NextKeyline(fl,WhiteSpace); //read title line and skip it
      line = _NextKeyline(fl,WhiteSpace); //read next directive line
      boost::to_upper(line);

      if( line.substr(0,4)=="UNIT" ) { //skip 'unit' line
	line = _NextKeyline(fl,WhiteSpace); //read next directive line
	boost::to_upper(line); 
      }

      if( line.substr(0,4)=="NEUT" ) { //skip 'neutral groups' line (DL_POLY Classic FIELD format)
	line = _NextKeyline(fl,WhiteSpace); //look for next directive line
	boost::to_upper(line);
      }

      int nmol_types;

      if( !_isKeyInt(line,WhiteSpace,"MOLEC",nmol_types) ) 
	throw std::runtime_error("Error: missing integer number in directive '" + line + "' in topology file '"+ _fname +"'");

#ifdef DEBUG
      cout << "Read from topology file " << _fname << " : '" << line << "' - " << nmol_types << endl;
#endif

      string mol_name;

      int id=0;
      for (int nmol_type=0;nmol_type<nmol_types; nmol_type++) {

	mol_name = _NextKeyline(fl,WhiteSpace);
        Molecule *mi = top.CreateMolecule(mol_name);

	int nreplica = 1;
	line = _NextKeyInt(fl,WhiteSpace,"NUMMOL",nreplica);

#ifdef DEBUG
	cout << "Read from topology file " << _fname << " : '" << mol_name << "' - '" << line << "' - " << nreplica << endl;
#endif

	line = _NextKeyInt(fl,WhiteSpace,"ATOMS",natoms);

#ifdef DEBUG
	cout << "Read from topology file " << _fname << " : '" << line << "' - " << natoms << endl;
#endif

	//read molecule
	int id_map[natoms];
	for (int i=0;i<natoms;){//i is altered in repeater loop
	  string beadtype;
	  fl >> beadtype;
	  BeadType *type = top.GetOrCreateBeadType(beadtype);
	  double mass;
	  fl >> mass;
	  double charge;
	  fl >> charge;

          getline(fl,line); //rest of the atom line

          Tokenizer tok(line, WhiteSpace);
	  vector<string> fields;
	  tok.ToVector(fields);

	  int repeater=1;
	  if (fields.size() > 1) repeater=boost::lexical_cast<int>(fields[0]);

	  for(int j=0;j<repeater;j++){

	    string beadname = beadtype + "#" + boost::lexical_cast<string>(i+1);
	    Bead *bead = top.CreateBead(1, beadname, type, res->getId(), mass, charge);

            stringstream nm;
            nm << bead->getResnr() + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
            mi->AddBead(bead, nm.str());
	    id_map[i]=bead->getId();
	    i++;
#ifdef DEBUG
	    cout << "Atom identification in maps '" << nm.str() << "'" << endl;
#endif
	  }
	  matoms += repeater;
	}

	while (line != "FINISH"){
	  if ((line == "BONDS")||(line == "ANGLES")||(line == "DIHEDRALS")) {
	    string type = line;
	    int count;
	    fl >> count;
	    for (int i=0;i<count;i++){
	      fl >> line; //bond/angle/dih type not used
	      int ids[4];
              Interaction *ic;
	      fl >> ids[0]; fl>>ids[1];
	      if (type == "BONDS"){
	        ic = new IBond(id_map[ids[0]-1],id_map[ids[1]-1]); // -1 due to fortran vs c 
	      } else if (type == "ANGLES"){
		fl >> ids[2];
	        ic = new IAngle(id_map[ids[0]-1],id_map[ids[1]-1],id_map[ids[2]-1]); // -1 due to fortran vs c 
	      } else if (type == "DIHEDRALS"){
		fl >> ids[2]; fl >> ids[3];
	        ic = new IDihedral(id_map[ids[0]-1],id_map[ids[1]-1],id_map[ids[2]-1],id_map[ids[3]-1]); // -1 due to fortran vs c 
	      }
              ic->setGroup(type);
              ic->setIndex(i);
              ic->setMolecule(mi->getId());
              top.AddBondedInteraction(ic);
              mi->AddInteraction(ic);
	      getline(fl,line);
	    }
	  }
	  fl >> line;
	  boost::to_upper(line);
	  if (fl.eof())
            throw std::runtime_error("Error: unexpected end of dlpoly file " + _fname + " while scanning for kerword 'finish'");
	}

#ifdef DEBUG
	cout << "Read from topology file " << _fname << " : '" << line << "' - done with '" << mol_name << "'" << endl;
#endif

	getline(fl, line); //rest of the FINISH line

	//replicate molecule
	for (int replica=1;replica<nreplica;replica++){
          Molecule *mi_replica = top.CreateMolecule(mol_name);
	  for(int i=0;i<mi->BeadCount();i++){
	    Bead *bead=mi->getBead(i);
	    BeadType *type = top.GetOrCreateBeadType(bead->Type()->getName());
	    string beadname=mi->getBeadName(i);
	    Bead *bead_replica = top.CreateBead(1, bead->getName(), type, res->getId(), bead->getM(), bead->getQ());
	    mi_replica->AddBead(bead_replica,beadname);
	  }
	  matoms+=mi->BeadCount();
	  InteractionContainer ics=mi->Interactions();
          for(vector<Interaction *>::iterator ic=ics.begin(); ic!=ics.end(); ++ic) {
            Interaction *ic_replica;
	    int offset = mi_replica->getBead(0)->getId() - mi->getBead(0)->getId();
	    if ((*ic)->BeadCount() == 2) {
	      ic_replica = new IBond((*ic)->getBeadId(0)+offset,(*ic)->getBeadId(1)+offset);
	    } else if ((*ic)->BeadCount() == 3) {
	      ic_replica = new IAngle((*ic)->getBeadId(0)+offset,(*ic)->getBeadId(1)+offset,(*ic)->getBeadId(2)+offset);
	    } else if ((*ic)->BeadCount() == 4) {
	      ic_replica = new IDihedral((*ic)->getBeadId(0)+offset,(*ic)->getBeadId(1)+offset,(*ic)->getBeadId(2)+offset,(*ic)->getBeadId(3)+offset);
	    }
            ic_replica->setGroup((*ic)->getGroup());
            ic_replica->setIndex((*ic)->getIndex());
            ic_replica->setMolecule(mi_replica->getId());
            top.AddBondedInteraction(ic_replica);
            mi_replica->AddInteraction(ic_replica);
	  }
	}
      }
      top.RebuildExclusions();
    }

#ifdef DEBUG
    getline(fl, line); //is "close" found?
    if(line=="close") {
      cout << "Read from topology file " << _fname << " : '" << line << "' - done with topology" << endl;
    }
    else {
      	    cout << "Read from topology file " << _fname << " : 'EOF' - done with topology (directive 'close' not read!)" << endl;
    }
#endif

    //we don't need the rest
    fl.close();

    string filename;

    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        filename="CONFIG";
      } else {
	filename=filepath.parent_path().string() + "/CONFIG";
      }
    } else {
      filename=filepath.parent_path().string()+boost::filesystem::basename(filepath)+".dlpc";
      cout << "NOTE: explicit dlpoly topology file name given '" << file << "', so trying to read boundary conditions from CONFIG file named '" << filename << "'" << endl;
    }

    fl.open(filename.c_str());

    if(fl.is_open()) {
      string line;

      getline(fl, line); //title

#ifdef DEBUG
      cout << "Read from initial configuration file : '" << line << "' - header" << endl;
#endif

      getline(fl, line); // 2nd header line

#ifdef DEBUG
      cout << "Read from initial configuration file : '" << line << "' - directives line" << endl;
#endif

      Tokenizer tok(line, WhiteSpace);
      vector<string> fields;
      tok.ToVector(fields);

      if( fields.size() < 3 ) 
	throw std::runtime_error("Error: too few directive switches (<3) in the initial configuration (check its 2-nd line)");	

      mavecs = boost::lexical_cast<int>(fields[0]);
      mpbct  = boost::lexical_cast<int>(fields[1]);
      natoms = boost::lexical_cast<int>(fields[2]);

      hasVs = (mavecs > 0); // 1 or 2 => in DL_POLY frame velocity vector follows coords for each atom/bead
      hasFs = (mavecs > 1); // 2      => in DL_POLY frame force vector follows velocities for each atom/bead

      top.SetHasVel(hasVs);
      top.SetHasForce(hasFs);

      if(natoms != matoms)
#ifdef DEBUG
	cout << "Warning: N of atoms/beads in initial configuration & topology differ: " << natoms << " =?= " << matoms << endl;
#else
      	throw std::runtime_error("Error: N of atoms/beads in initial configuration & topology differ " +
	    boost::lexical_cast<string>(natoms) + " vs " + boost::lexical_cast<string>(matoms));
#endif

      vec box_vectors[3];
      for (int i=0;i<3;i++){ // read 3 box lines
        getline(fl, line);

#ifdef DEBUG
	cout << "Read from initial configuration file : '" << line << "' - box vector # " << i+1 << " (Angs)" << endl;
#endif

        if(fl.eof())
          throw std::runtime_error("Error: unexpected EOF in dlpoly file " + filename +", when reading box vector " +
              boost::lexical_cast<string>(i));

        Tokenizer tok(line, WhiteSpace);
        vector<double> fields;
        tok.ConvertToVector<double>(fields);
        box_vectors[i]=vec(fields[0],fields[1],fields[2]);

#ifdef DEBUG
	cout << "Read from initial configuration file : '" << fixed << setprecision(10) << setw(20) << fields[0] << setw(20) << fields[1] << setw(20) << fields[2] << "' - box vector # " << i+1 << " (Angs)" << endl;
#endif
      }
      matrix box(box_vectors[0],box_vectors[1],box_vectors[2]);

      if(mpbct == 0) {
	pbc_type=BoundaryCondition::typeOpen;
      }
      else if(mpbct == 1 || mpbct == 2 ) {
	pbc_type=BoundaryCondition::typeOrthorhombic;
      }
      else if(mpbct == 3) {
	pbc_type=BoundaryCondition::typeTriclinic;
      }

      top.setBox(box,pbc_type);

      fl.close();

#ifdef DEBUG
      cout << "Read from initial configuration file : box/cell matrix - done with boundaries" << endl << endl;
#endif
    }
    else {
      cout << "NOTE: could not open dlpoly file " << filename << ", so no PBC set in topology - assuming 'open box'" << endl;
    }

    return true;
}

}}

