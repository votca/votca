/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/lexical_cast.hpp>
#include <votca/tools/getline.h>
#include <boost/algorithm/string.hpp>
#include "grotopologyreader.h"


namespace votca { namespace csg {

bool GROTopologyReader::ReadTopology(string file, Topology &top)
{ 
    // cleanup topology to store new data
    ifstream fl;
    string tmp;
    top.Cleanup();
    
    fl.open(file.c_str());
    if(!fl.is_open())
      throw std::runtime_error("Could not open file: " + file);

    string title;
    getline(fl, title);
    
    getline(fl, tmp);
    int natoms = atoi(tmp.c_str());
    for(;natoms>0; natoms--) {
        char c[6];
        fl.read(c, 5);
        c[5] = 0; //fix trailing char

        int resnr = atoi(c);
        //residue name
        fl.read(c, 5);
        if(resnr >= top.ResidueCount()) {
            if (top.ResidueCount()==0) //gro resnr start with 1 but VOTCA starts with 0
              top.CreateResidue("ZERO"); // create 0 res, to allow to keep gro numbering
            string withoutWhitespace;
            withoutWhitespace  = string(c);
            boost::algorithm::trim(withoutWhitespace);
            top.CreateResidue(withoutWhitespace);
        }
        //atom name
        char atomname[6];
        fl.read(atomname, 5);
        atomname[5]=0; //fix trailing character
        string atomnameWoWs;
        atomnameWoWs  = string(atomname);
        boost::algorithm::trim(atomnameWoWs);
        //atom number, not needed
        char buf[5];
        fl.read(buf, 5);
        //position
        char x[8], y[8], z[8];
        fl.read(x,8);
        fl.read(y,8);
        fl.read(z,8);
        top.CreateBead(1, atomnameWoWs, top.GetOrCreateBeadType(atomnameWoWs), resnr, 1., 0.);
        getline(fl, tmp);
    }

    getline(fl, tmp); //read box line
    if(fl.eof())
      throw std::runtime_error("unexpected end of file in poly file, when boxline");
    Tokenizer tok(tmp, " ");
    vector<double> fields;
    tok.ConvertToVector<double>(fields);
    matrix box;
    if ( fields.size() == 3 ) {
      box.ZeroMatrix();
      for (int i=0;i<3;i++) box[i][i] = fields[i];
    } else if ( fields.size() == 9 ) {
      box[0][0]=  fields[0];
      box[1][1]=  fields[1];
      box[2][2]=  fields[2];
      box[0][1]=  fields[3];
      box[0][2]=  fields[4];
      box[1][1]=  fields[5];
      box[1][2]=  fields[6];
      box[2][1]=  fields[7];
      box[2][2]=  fields[8];
    } else {
      throw std::runtime_error("Error while reading box (last) line of"+ file);
    }
    top.setBox(box);

    fl.close();
    
    cout << "WARNING: topology created from .gro file, masses and charges are wrong!\n";
    
    return true;
}

}}
