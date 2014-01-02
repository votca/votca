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

#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/tools/getline.h>
#include "pdbreader.h"

namespace votca { namespace csg {
    using namespace boost;
using namespace std;

bool PDBReader::ReadTopology(string file,  Topology &top)
{
   _topology = true;
   top.Cleanup();

   _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topologyl file: " + file);

   NextFrame(top);

    _fl.close();

    return true;
}

bool PDBReader::Open(const string &file)
{
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topologyl file: " + file);
    return true;
}

void PDBReader::Close()
{
    _fl.close();
}

bool PDBReader::FirstFrame(Topology &top)
{
    _topology = false;
    NextFrame(top);
    return true;
}

bool PDBReader::NextFrame(Topology &top)
{
    string line;
    int i = 0 ;
    while ( std::getline(_fl, line) ){
        if( wildcmp("ATOM*",line.c_str()) || wildcmp("HETATM*",line.c_str())){
            
            //      according to PDB format
	    string x,y,z, resNum, resName, atName;
            try {
	      /* Some pdb don't include all this, read only what we realy need*/
	      /* leave this here in case we need more later*/
              //string recType    (line,( 1-1),6); // str,  "ATOM", "HETATM"
              //string atNum      (line,( 7-1),6); // int,  Atom serial number
              atName=string(line,(13-1),4); // str,  Atom name
              //string atAltLoc   (line,(17-1),1); // char, Alternate location indicator
              resName=string(line,(18-1),4); // str,  Residue name
              //string chainID    (line,(22-1),1); // char, Chain identifier
              resNum=string(line,(23-1),4); // int,  Residue sequence number
              //string atICode    (line,(27-1),1); // char, Code for insertion of res
              x=string(line,(31-1),8); // float 8.3 ,x
              y=string(line,(39-1),8); // float 8.3 ,y
              z=string(line,(47-1),8); // float 8.3 ,z
              //string atOccup    (line,(55-1),6); // float  6.2, Occupancy
              //string atTFactor  (line,(61-1),6); // float  6.2, Temperature factor
              //string segID      (line,(73-1),4); // str, Segment identifier
              //string atElement  (line,(77-1),2); // str, Element symbol
              //string atCharge   (line,(79-1),2); // str, Charge on the atom

	    } catch (std::out_of_range& err) {
	      throw std::runtime_error("Misformated pdb file");
	    }
            boost::algorithm::trim(atName);
            boost::algorithm::trim(resName);
            boost::algorithm::trim(resNum);
            boost::algorithm::trim(x);
            boost::algorithm::trim(y);
            boost::algorithm::trim(z);

	    i++;
            if(!_topology && i > top.BeadCount())
              throw std::runtime_error("number of beads in topology and trajectory differ");

            Bead *b;
            if(_topology){
	      int resnr = boost::lexical_cast<int>(resNum);
              if(resnr >= top.ResidueCount()) {
                if (top.ResidueCount()==0) //pdb resnr start with 1 but VOTCA starts with 0
	          top.CreateResidue("ZERO"); // create 0 res, to allow to keep pdb numbering
                top.CreateResidue(resName);
	      }
              //this is not correct, but still better than no type at all!
	      BeadType *type = top.GetOrCreateBeadType(atName);
              
	      b = top.CreateBead(1, atName, type, resnr, 1., 0.);
	    } else {
                b = top.getBead(i-1);
	    }

            b->setPos(vec(
                    boost::lexical_cast<double>(x),
                    boost::lexical_cast<double>(y),
                    boost::lexical_cast<double>(z)
                ));

	}

        if (( line == "ENDMDL" )||( _fl.eof())){
	  break;
	}
    }

    if(!_topology && (i>0) && i != top.BeadCount())
      throw std::runtime_error("number of beads in topology and trajectory differ");
   
    if (_topology)
      cout << "WARNING: topology created from .pdb file, charges and masses are wrong!\n";
    
    return !_fl.eof();
}

}}

