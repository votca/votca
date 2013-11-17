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
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp> 
#include <votca/tools/getline.h>

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "dlpolytopologyreader.h"
#ifdef DLPOLY_FORTRAN
#include "dlpoly/dlp_io_layer.h"
#include "fortan_mangling.h"
#endif

namespace votca { namespace csg {

bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top)
{
    std::ifstream fl;
    boost::filesystem::path filepath(file.c_str());
    string filename;
#ifdef DLPOLY_FORTRAN
    if (file != ".dlpoly")
      throw std::runtime_error("Reading from different filename/directories not implemented yet. (use --top '.dlpoly')");

    struct FieldSpecsT  FieldBase;
    struct FrameSpecsT  FrameBase;
    struct MolecSpecsT *MolecBase;
    struct FieldSiteT  *FieldSite;
    struct FrameSiteT  *FrameSite;

    int idnode,matms,natms,nmols,nmolt;
    int istateF;

    int inode=matms=natms=nmols=nmolt=0;

    // TODO: istateF must be an enum!
    istateF=1;

    // TODO: we need to fix the file naming!
    FortranCInterface_GLOBAL(field_scan,FIELD_SCAN)(&istateF,&matms,&natms,&nmolt);

    MolecBase = new MolecSpecsT[nmolt];
    FieldSite = new FieldSiteT[natms];

    FieldBase.nmols = nmolt;
    FieldBase.natms = natms;

    FortranCInterface_GLOBAL(field_read,FIELD_READ)(&istateF,&FieldBase,MolecBase,FieldSite);

    // AB: if on return istateF < 0  => in the next F-call the relevant F-arrays will be deallocated (at the end)
    // AB: NOT TO RE-/DE-ALLOCATE F-arrays in the next F-call, reset istateF = 0
    istateF = 0;

    // TODO: fix residue naming / assignment
    Residue *res = top.CreateResidue("no");

    // read the atoms
    int mol_offset=0;
    for(int im=0; im<nmolt; im++){
        for(int imr=0; imr<MolecBase[im].nrept; ++imr) {
            Molecule *mi = top.CreateMolecule(MolecBase[im].name);
            for(int ims=0; ims<MolecBase[im].nsites; ims++) {
	        int is=mol_offset+ims;
                BeadType *type = top.GetOrCreateBeadType(FieldSite[is].type); // what is
	        string beadname = boost::lexical_cast<string>(FieldSite[is].name) + "#" + boost::lexical_cast<string>(ims+1);
                Bead *bead = top.CreateBead(1, beadname, type, res->getId(), FieldSite[is].m, FieldSite[is].q);

                stringstream nm;
                nm << bead->getResnr() + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
                mi->AddBead(bead, nm.str());
            }
        }
	mol_offset+=MolecBase[im].nsites;
    }

    delete [] MolecBase;
    delete [] FieldSite;
#else
    // TODO: fix residue naming / assignment
    Residue *res = top.CreateResidue("no");

    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        filename="FIELD";
      } else {
	filename=filepath.parent_path().string() + "/FIELD";
      }
    } else {
      filename=file;
    }
    fl.open(filename.c_str());
    if (!fl.is_open()){
      throw std::runtime_error("could open dlpoly file " + filename);
    } else {
      string line;
      getline(fl, line); //title
      getline(fl, line); //unit line
      fl >> line;
      if (line != "MOLECULES")
          throw std::runtime_error("unexpected line in dlpoly file " + filename + ", expected 'MOLECULES' got" + line);
      int nmol_types;
      fl >> nmol_types;
      getline(fl, line); //rest of MOLECULES line
      int id=0;
      for (int nmol_type=0;nmol_type<nmol_types; nmol_type++){
	string mol_name;
        getline(fl, mol_name); //molecule name might incl. spaces
	boost::erase_all(mol_name, " ");
        Molecule *mi = top.CreateMolecule(mol_name);
        fl >> line;
        if (line != "NUMMOLS")
          throw std::runtime_error("unexpected line in dlpoly file " + filename + ", expected 'NUMMOLS' got" + line);
	int nreplica;
	fl >> nreplica;
	fl >> line;
        if (line != "ATOMS")
          throw std::runtime_error("unexpected line in dlpoly file " + filename + ", expected 'ATOMS' got" + line);
	int natoms;
	fl >> natoms;
	//read molecule
	int id_map[natoms];
	for (int i=0;i<natoms;){//i is altered in reapeater loop
	  string beadtype;
	  fl >> beadtype;
	  BeadType *type = top.GetOrCreateBeadType(beadtype);
	  double mass;
	  fl >> mass;
	  double charge;
	  fl >> charge;
          getline(fl,line); //rest of the atom line
          Tokenizer tok(line, " ");
          vector<int> fields;
          tok.ConvertToVector<int>(fields);
	  int repeater=1;
	  if (fields.size() > 1)
	    repeater=fields[0];
	  for(int j=0;j<repeater;j++){
	    string beadname = beadtype + "#" + boost::lexical_cast<string>(j+1);
	    Bead *bead = top.CreateBead(1, beadname , type, res->getId(), mass, charge);
            stringstream nm;
            nm << bead->getResnr() + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
            mi->AddBead(bead, nm.str());
	    id_map[i]=bead->getId();
	    i++;
	  }
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
	  if (fl.eof())
            throw std::runtime_error("unexpected end of dlpoly file " + filename + " while scanning for 'FINISH'");
	}
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
	    //TODO copy interactions
	  }
	}
      }
    }
    //we don't need the rest.
    fl.close();
#endif
    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        filename="CONFIG";
      } else {
	filename=filepath.parent_path().string() + "/CONFIG";
      }
    } else {
      cout << "NOTE: Explicit dlploy topology filename given, so no config will be openend and no boundary conditions will set in the topology." << endl;
      return true;
    }
    fl.open(filename.c_str());
    if(fl.is_open()) {
      string line;
      getline(fl, line); //title
      getline(fl, line); // 2nd header line
      vec box_vectors[3];
      for (int i=0;i<3;i++){ // read 3 box lines
        getline(fl, line);
        if(fl.eof())
          throw std::runtime_error("unexpected end of file in dlpoly file " + filename +", when reading box vector"  +
              boost::lexical_cast<string>(i));
        Tokenizer tok(line, " ");
        vector<double> fields;
        tok.ConvertToVector<double>(fields);
        box_vectors[i]=vec(fields[0],fields[1],fields[2]);
      }
      matrix box(box_vectors[0],box_vectors[1],box_vectors[2]);
      top.setBox(box);
      fl.close();
    }
    else {
      cout << "NOTE: Could not open dlploy file " << filename << "so no boundary conditions, where set in the topology." << endl;
    }
    return true;
}

}}

