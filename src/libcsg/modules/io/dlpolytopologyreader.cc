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

#include "dlpolytopologyreader.h"
#include "dlpoly/dlp_io_layer.h"

namespace votca { namespace csg {

bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top)
{
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
    int is=0;
    for(int im=0; im<nmolt; im++){
        for(int imr=0; imr<MolecBase[im].nrept; ++imr) {
            Molecule *mi = top.CreateMolecule(MolecBase[im].name);
            for(int ims=0; ims<MolecBase[im].nsites; ims++) {

                BeadType *type = top.GetOrCreateBeadType(FieldSite[is].name); // what is
                Bead *bead = top.CreateBead(1, FieldSite[is].type, type, res->getId(), FieldSite[is].m, FieldSite[is].q);

                stringstream nm;
                nm << bead->getResnr() + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
                mi->AddBead(bead, nm.str());
            }
        }
    }

    delete [] MolecBase;
    delete [] FieldSite;

    std::ifstream fl;
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

    return true;
}

}}

