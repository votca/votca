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
#include <boost/filesystem/convenience.hpp>
#include "pdbtrajectoryreader.h"

namespace votca { namespace csg {
    using namespace boost;
using namespace std;

bool PDBTrajectoryReader::Open(const string &file)
{

    boost::filesystem::path filepath(file.c_str());
    
    // Check that the file has the correct extension
    if(boost::filesystem::extension(filepath).size()==0 ){
        throw std::ios_base::failure("ERROR on opening .pdb file'"
                     +file+"' - extension is expected, use .pdb");
    }else if(boost::filesystem::extension(filepath)!=".pdb"){
        throw std::ios_base::failure("Error on opening .pdb file '"
                     +file+"' - wrong extension, use .pdb");
    }

    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: " + file);
    return true;
}

void PDBTrajectoryReader::Close()
{
    _fl.close();
}

bool PDBTrajectoryReader::FirstFrame(Topology &top)
{
    _first_frame=true;
    bool res=NextFrame(top);
    _first_frame=false;
    return res;
}

bool PDBTrajectoryReader::NextFrame(Topology &top)
{
    string line;

    ////////////////////////////////////////////////////////////////////////////////
    // Read in information from .pdb file
    ////////////////////////////////////////////////////////////////////////////////
    int bead_count = 0 ;
    while ( std::getline(_fl, line) ){
        if( wildcmp("CRYST1*",line.c_str())){
             string a, b, c, alpha, beta, gamma;
             try {
                                            // 1 -  6       Record name    "CRYST1"
               a    =string(line,(7-1),9);  // 7 - 15       Real(9.3)      a           (Angstroms)
               b    =string(line,(16-1),9); //16 - 24       Real(9.3)      b           (Angstroms)
               c    =string(line,(25-1),9); //25 - 33       Real(9.3)      c           (Angstroms)
               alpha=string(line,(34-1),7); //34 - 40       Real(7.2)      alpha       (degrees)
               beta =string(line,(41-1),7); //41 - 47       Real(7.2)      beta        (degrees)
               gamma=string(line,(48-1),7); //48 - 54       Real(7.2)      gamma       (degrees)
                                            //56 - 66       LString        Space group
                                            //67 - 70       Integer        Z value
            } catch (std::out_of_range& err) {
              throw std::runtime_error("Misformated pdb file in CRYST1 line");
            }
            boost::algorithm::trim(a);
            boost::algorithm::trim(b);
            boost::algorithm::trim(c);
            boost::algorithm::trim(alpha);
            boost::algorithm::trim(beta);
            boost::algorithm::trim(gamma);
	        if ((!wildcmp("90*",alpha.c_str()))||(!wildcmp("90*",alpha.c_str()))||(!wildcmp("90*",alpha.c_str()))){
	         throw std::runtime_error("Non cubical box in pdb file not implemented, yet!");
            }
            double aa = boost::lexical_cast<double>(a)/10.0;
            double bb = boost::lexical_cast<double>(b)/10.0;
            double cc = boost::lexical_cast<double>(c)/10.0;
	        top.setBox(matrix(vec(aa, 0 , 0 ),
	                          vec(0 , bb, 0 ),
                              vec(0 , 0 , cc)));

	    }
        if( wildcmp("ATOM*",line.c_str()) || wildcmp("HETATM*",line.c_str())){
            
            //      according to PDB format
	        string x,y,z, resNum, resName, atName;
            try {
	      /* Some pdb don't include all this, read only what we realy need*/
	      /* leave this here in case we need more later*/
              //string recType    (line,( 1-1),6); // str       ,  "ATOM", "HETATM"
              //string atNum      (line,( 7-1),6); // int       , Atom serial number
              atName =      string(line,(13-1),4); // str       , Atom name
              //string atAltLoc   (line,(17-1),1); // char      , Alternate location indicator
              resName=      string(line,(18-1),3); // str       , Residue name
              //string chainID    (line,(22-1),1); // char      , Chain identifier
              resNum=       string(line,(23-1),4); // int       , Residue sequence number
              //string atICode    (line,(27-1),1); // char      , Code for insertion of res
              x     =       string(line,(31-1),8); // float 8.3 , x
              y     =       string(line,(39-1),8); // float 8.3 , y
              z     =       string(line,(47-1),8); // float 8.3 , z
              //string atOccup    (line,(55-1),6); // float 6.2 , Occupancy
              //string atTFactor  (line,(61-1),6); // float 6.2 , Temperature factor
              //string segID      (line,(73-1),4); // str       , Segment identifier
              //string atElement  (line,(77-1),2); // str       , Element symbol
              //string atCharge   (line,(79-1),2); // str       , Charge on the atom

            } catch (std::out_of_range& err) {
              throw std::runtime_error("Misformated pdb file in atom line # "+ boost::lexical_cast<string>(bead_count));
            }
            boost::algorithm::trim(atName);
            boost::algorithm::trim(resName);
            boost::algorithm::trim(resNum);
            boost::algorithm::trim(x);
            boost::algorithm::trim(y);
            boost::algorithm::trim(z);

            bead_count++;
            if(bead_count > top.BeadCount())
              throw std::runtime_error("number of beads in topology and trajectory differ");

            Bead *b;
            b = top.getBead(bead_count-1);
            // convert to nm from A
            b->setPos(vec(
                        boost::lexical_cast<double>(x)/10.0,
                        boost::lexical_cast<double>(y)/10.0,
                        boost::lexical_cast<double>(z)/10.0
                        ));

        }

        if (( line == "ENDMDL" ) || ( line == "END" ) || ( _fl.eof())){
            break;
        }
    }

    if((bead_count>0) && bead_count != top.BeadCount())
      throw std::runtime_error("number of beads in topology and trajectory differ");
   
    return !_fl.eof();
}

}}

