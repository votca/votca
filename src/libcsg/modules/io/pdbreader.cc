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

   if(_topology)
        top.CreateResidue("DUM");

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
    int natoms = 0 ;
    while ( std::getline(_fl, line,'\n') ){
        if(     boost::find_first(line, "ATOM"  )   || 
                boost::find_first(line, "HETATM")     //TODO I think this should be wildcmp ATOM*
		//other ATOM in the middle of the line will trigger this, too.
                ){
            
            //      according to PDB format
            string recType    (line,( 1-1),6); // str,  "ATOM", "HETATM"
            string atNum      (line,( 7-1),6); // int,  Atom serial number
            string atName     (line,(13-1),4); // str,  Atom name
            string atAltLoc   (line,(17-1),1); // char, Alternate location indicator
            string resName    (line,(18-1),4); // str,  Residue name
            string chainID    (line,(22-1),1); // char, Chain identifier
            string resNum     (line,(23-1),4); // int,  Residue sequence number
            string atICode    (line,(27-1),1); // char, Code for insertion of res
            string x          (line,(31-1),8); // float 8.3 ,x
            string y          (line,(39-1),8); // float 8.3 ,y
            string z          (line,(47-1),8); // float 8.3 ,z
            string atOccup    (line,(55-1),6); // float  6.2, Occupancy
            string atTFactor  (line,(61-1),6); // float  6.2, Temperature factor
            string segID      (line,(73-1),4); // str, Segment identifier
            string atElement  (line,(77-1),2); // str, Element symbol
            string atCharge   (line,(79-1),2); // str, Charge on the atom

	    boost::algorithm::trim(recType);
            algorithm::trim(atNum);
            algorithm::trim(atName);
            algorithm::trim(atAltLoc);
            algorithm::trim(resName);
            algorithm::trim(chainID);
            algorithm::trim(resNum);
            algorithm::trim(atICode);
            algorithm::trim(x);
            algorithm::trim(y);
            algorithm::trim(z);
            algorithm::trim(atOccup);
            algorithm::trim(atTFactor);
            algorithm::trim(segID);
            algorithm::trim(atElement);
            algorithm::trim(atCharge);

	    natoms++;
            Bead *b;
            if(_topology){
	         //TODO
                //b = top.CreateBead(1, fields[0]+boost::lexical_cast<string>(i),
                //        top.GetOrCreateBeadType(fields[0]), 0, 0, 0);
	    } else {
                b = top.getBead(natoms);
	    }

            b->setPos(vec(
                    boost::lexical_cast<double>(x),
                    boost::lexical_cast<double>(y),
                    boost::lexical_cast<double>(z)
                ));

	}
    }
    
    if(!_topology && natoms !=top.BeadCount())
      throw std::runtime_error("number of beads in topology and trajectory differ");

    return !_fl.eof();;
}

}}

