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

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "dlpolytopologyreader.h"
#ifdef DLPOLY
#include "dlpoly/dlp_io_layer.h"
#endif

namespace votca { namespace csg {

bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top)
{
    std::ifstream fl;
#ifdef DLPOLY
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
    field_scan_(&istateF,&matms,&natms,&nmolt);

    MolecBase = new MolecSpecsT[nmolt];
    FieldSite = new FieldSiteT[natms];

    FieldBase.nmols = nmolt;
    FieldBase.natms = natms;

    field_read_(&istateF,&FieldBase,MolecBase,FieldSite);

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

    fl.open("FIELD");
    if (!fl.is_open()){
      throw std::runtime_error("could open dlpoly file FIELD");
    } else {
      string line;
      getline(fl, line); //title
      getline(fl, line); //unit line
      fl >> line;
      if (line != "MOLECULES")
          throw std::runtime_error("unexpected line in dlpoly file FIELD, expected 'MOLECULES' got" + line);
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
          throw std::runtime_error("unexpected line in dlpoly file FIELD, expected 'NUMMOLS' got" + line);
	int nreplica;
	fl >> nreplica;
	fl >> line;
        if (line != "ATOMS")
          throw std::runtime_error("unexpected line in dlpoly file FIELD, expected 'ATOMS' got" + line);
	int natoms;
	fl >> natoms;
	//read molecule
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
	    i++;
	  }
	}
	while (line != "FINISH"){
	  fl >> line;
	  if (fl.eof())
            throw std::runtime_error("unexpected end of dlpoly file FIELD while scanning for 'FINISH'");
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
	  }
	}
      }
    }
    //we don't need the rest.
    fl.close();
#endif
    fl.open("CONFIG");
    if(fl.is_open()) {
      string line;
      getline(fl, line); //title
      getline(fl, line); // 2nd header line
      vec box_vectors[3];
      for (int i=0;i<3;i++){ // read 3 box lines
        getline(fl, line);
        if(fl.eof())
          throw std::runtime_error("unexpected end of file in dlpoly file CONFIG, when reading box vector"  +
              boost::lexical_cast<string>(i));
        Tokenizer tok(line, " ");
        vector<double> fields;
        tok.ConvertToVector<double>(fields);
        box_vectors[i]=vec(fields[0],fields[1],fields[2]);
      }
      matrix box(box_vectors[0],box_vectors[1],box_vectors[2]);
      top.setBox(box);
    }
    else {
      cout << "NOTE: Could open dlploy file CONFIG, so no boundary conditions, where set in the topology" << endl;
    }
    fl.close();

    return true;
}

}}

