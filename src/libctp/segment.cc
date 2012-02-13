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
Segment::Segment(int id, string name)
        : _id(id), _name(name), _hasOccProb(0),
          _hasLambdas(0) { }

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

    _lambdasIntra.clear();
    _eMpoles.clear();
    _occProb.clear();

}


void Segment::setOcc(int carrier, double occ) {
    _hasOccProb = true;
    _occProb[carrier] = occ;
}
const double &Segment::getOcc(int carrier) {
    return _occProb.at(carrier);
}

void Segment::setESiteIntra(int carrier, double energy) {
    _eSiteIntra[carrier] = energy;
    _hasESiteIntra = true;
}
const double &Segment::getESiteIntra(int carrier) {
    return _eSiteIntra[carrier];
}

void Segment::setLambdaIntra(int state0, int state1, double lambda) {
    _hasLambdas = true;
    _lambdasIntra[state0][state1] = lambda;
}
const double &Segment::getLambdaIntra(int state0, int state1) {
    return _lambdasIntra.at(state0).at(state1); 
}

void Segment::setEMpoles(int state, double energy) {
    _hasChrgState[state] = true;
    _eMpoles[state] = energy;
}
const double &Segment::getEMpoles(int state) {
    return _eMpoles.at(state); 
}


void Segment::AddChrgState(int state, bool yesno) {
    this->_hasChrgState[state] = yesno;
}

void Segment::chrg(int state) {
    vector < Atom* > ::iterator ait;
    for (ait = this->Atoms().begin();
            ait < this->Atoms().end();
            ait++) {
        (*ait)->chrg(state);
    }
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




void Segment::Rigidify() {

    if (this->getType()->canRigidify()) {
        // Establish which atoms to use to define local frame
        vector<Fragment*> ::iterator fit;

        for (fit = this->Fragments().begin();
                fit < this->Fragments().end();
                fit++) {
                (*fit)->Rigidify();
        }
    }
    else { return; }
}















void Segment::WritePDB(FILE *out, string tag1, string tag2) {

  if (tag1 == "Fragments") {
    vector < Fragment* > :: iterator frag;
    for (frag = _fragments.begin(); frag < _fragments.end(); ++frag){
         int id = (*frag)->getId();
         string name =  (*frag)->getName();
         name.resize(3);
         string resname = (*frag)->getSegment()->getName();
         resname.resize(3);
         int resnr = (*frag)->getSegment()->getId();
         vec position = (*frag)->getPos();  

         fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
                      "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n",
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
  if ( tag1 == "Atoms") {
    vector < Atom* > :: iterator atm;
    for (atm = _atoms.begin(); atm < _atoms.end(); ++atm) {
         int id = (*atm)->getId();
         string name =  (*atm)->getName();
         name.resize(3);
         string resname = (*atm)->getResname();
         resname.resize(3);
         int resnr = (*atm)->getResnr();
         vec position;
         if (tag2 == "MD")      { position = (*atm)->getPos(); }
         else if (tag2 == "QM") { position = (*atm)->getQMPos(); }

         fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
                      "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n",
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
}


}}
