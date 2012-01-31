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

#include <votca/ctp/segment.h>

using namespace std;

namespace votca { namespace ctp {
   
/// Default constructor
Segment::Segment(int id, string name) { 
    _id = id; _name = name; _occ = -1;
}

/// Destructor
Segment::~Segment() {

    vector < Fragment* > ::iterator fragit;
    for (fragit = this->Fragments().begin();
            fragit < this->Fragments().end();
            fragit++) {
        delete *fragit;
    }
    _fragments.clear();
    _atoms.clear();
}

/// Returns the ID of the molecule
const int &Segment::getId() {
    return _id;
}

/// Returns the name of the molecule
const string &Segment::getName() {
    return _name;
}

void Segment::AddFragment( Fragment* fragment ) {
    _fragments.push_back( fragment );
    fragment->setSegment(this);
}

void Segment::AddAtom( Atom* atom ) {
    _atoms.push_back( atom );
    atom->setSegment(this);
}


void Segment::calcPos() {
    vec pos = vec(0,0,0);
    double totWeight = 0.0;

    for (int i = 0; i< _atoms.size(); i++) {
        pos += _atoms[i]->getPos() * _atoms[i]->getWeight();
        totWeight += _atoms[i]->getWeight();
    }

    _CoM = pos / totWeight;
}


void Segment::WritePDB(FILE *out) {

    vector < Fragment* > :: iterator frag;
    for (frag = _fragments.begin(); frag < _fragments.end(); ++frag){
         int id = (*frag)->getId();
         string name =  (*frag)->getName();
         name.resize(3);
         string resname = (*frag)->getSegment()->getName();
         resname.resize(3);
         int resnr = (*frag)->getSegment()->getId();
         vec position = (*frag)->getPos();  

         fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n",
                 id,                    // Atom serial number           %5d
                 name.c_str(),          // Atom name                    %4s
                 " ",                   // alternate location indicator.%1s
                 resname.c_str(),       // Residue name.                %3s
                 "A",                   // Chain identifier             %1s
                 resnr,                 // Residue sequence number      %4d
                 " ",                   // Insertion of residues.       %1s
                 position.getX()*10,    // X in Angstroms               %8.3f
                 position.getY()*10,    // Y in Angstroms               %8.3f
                 position.getZ()*10,    // Z in Angstroms               %8.3f
                 1.0,                   // Occupancy                    %6.2f
                 0.0,                   // Temperature factor           %6.2f
                 " ",                   // Segment identifier           %4s
                 name.c_str(),          // Element symbol               %2s
                 " "                    // Charge on the atom.          %2s
                 );

    }

}


}}
