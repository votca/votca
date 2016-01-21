/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/segment.h>

using namespace std;

namespace votca { namespace xtp {
   
/// Default constructor
Segment::Segment(int id, string name)
        : _id(id),        _name(name),
          _has_e(false),  _has_h(false),_has_s(false),  _has_t(false)
            { _eMpoles.resize(5); }

// This constructor creates a copy of the stencil segment, without
// adding it to any containers further up in the hierarchy; i.e. the topology
// and molecules will neither know about the existence of this segment, nor
// be able to access it. Used for creating the ghost in PB corrected pairs.
Segment::Segment(Segment *stencil)
        : _id(stencil->getId()),    _name(stencil->getName()+"_ghost"),
          _typ(stencil->getType()), _top(NULL), _mol(NULL),
          _CoM(stencil->getPos()),
          _has_e(false), _has_h(false),_has_s(false),  _has_t(false)
            { _eMpoles.resize(5);

    vector<Fragment*> ::iterator fit;
    for (fit = stencil->Fragments().begin();
         fit < stencil->Fragments().end();
            fit++) {

        Fragment *newFrag = new Fragment(*fit);
        this->AddFragment(newFrag);

        vector<Atom*> ::iterator ait;
        for (ait = newFrag->Atoms().begin();
             ait < newFrag->Atoms().end();
             ait++) {
            this->AddAtom(*ait);
        }
    }
}


Segment::~Segment() {

    vector < Fragment* > ::iterator fragit;
    for (fragit = this->Fragments().begin();
            fragit < this->Fragments().end();
            fragit++) {
        delete *fragit;
    }
    _fragments.clear();
    _atoms.clear();

    _eMpoles.clear();
    _polarSites.clear();
}

void Segment::TranslateBy(const vec &shift) {

    _CoM = _CoM + shift;

    vector < Fragment* > ::iterator fit;
    for (fit = _fragments.begin();
            fit < _fragments.end();
            fit++) {

        (*fit)->TranslateBy(shift);
    }
}


void Segment::setHasState(bool yesno, int state) {

    if (state == -1) {
        _has_e = yesno;
    }
    else if (state == +1) {
        _has_h = yesno;
    }
    else if (state == +2) {
        _has_s = yesno;
    }
    else if (state == +3) {
        _has_t = yesno;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00e11h__");
    }
}

bool Segment::hasState(int state) {
    bool result;
    if (state == -1) {
        result = _has_e;
    }
    else if (state == +1) {
        result = _has_h;
    }
     else if (state == +2) {
        result = _has_s;
    }
     else if (state == +3) {
        result = _has_t;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00s11o__");
    }
    return result;
}

void Segment::setOcc(double occ, int e_h_s_t) {
    
    if (e_h_s_t == -1) {
        _occ_e = occ;
    }
    else if (e_h_s_t == +1) {
        _occ_h = occ;
    }
     else if (e_h_s_t == +2) {
        _occ_s = occ;
    }
     else if (e_h_s_t == +3) {
        _occ_t = occ;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00s11o__");
    }
}

    double Segment::getOcc(int e_h_s_t) {
    double result;
    if (e_h_s_t == -1) {
        result=_occ_e;
    }
    else if (e_h_s_t == +1) {
        result=_occ_h;
    }
     else if (e_h_s_t == +2) {
        result=_occ_s;
    }
     else if (e_h_s_t == +3) {
        result=_occ_t;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00s11o__"); // blabla what do I do here?
    }
    return result;
}
    
void Segment::setU_cC_nN(double dU, int state) {

    if (state == -1) {
        _U_cC_nN_e = dU;
    }
    else if (state == +1) {
        _U_cC_nN_h = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11a__");
    }
}


void Segment::setU_nC_nN(double dU, int state) {

    if (state == -1) {
        _U_nC_nN_e = dU;
    }
    else if (state == +1) {
        _U_nC_nN_h = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11b__");
    }
}


void Segment::setU_cN_cC(double dU, int state) {

    if (state == -1) {
        _U_cN_cC_e = dU;
    }
    else if (state == +1) {
        _U_cN_cC_h = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11c__");
    }
}

void Segment::setU_xX_nN(double dU, int state) {

    if (state == +2) {
        _U_xX_nN_s = dU;
    }
    else if (state == +3) {
        _U_xX_nN_t = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11d__"); //blabla?? What do I do here?
    }
}

void Segment::setU_nX_nN(double dU, int state) {

    if (state == +2) {
        _U_nX_nN_s = dU;
    }
    else if (state == +3) {
        _U_nX_nN_t = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11d__"); //blabla?? What do I do here?
    }
}

void Segment::setU_xN_xX(double dU, int state) {

    if (state == +2) {
        _U_xN_xX_s = dU;
    }
    else if (state == +3) {
        _U_xN_xX_t = dU;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00u11d__"); //blabla?? What do I do here?
    }
}

const double &Segment::getU_xX_nN(int state) {

    return (state == +3) ? _U_xX_nN_t : _U_xX_nN_s;
}

const double &Segment::getU_nX_nN(int state) {

    return (state == +3) ? _U_nX_nN_t : _U_nX_nN_s;
}

const double &Segment::getU_xN_xX(int state) {

    return (state == +3) ? _U_xN_xX_t : _U_xN_xX_s;
}



const double &Segment::getU_cC_nN(int state) {

    return (state == -1) ? _U_cC_nN_e : _U_cC_nN_h;
}


const double &Segment::getU_nC_nN(int state) {

    return (state == -1) ? _U_nC_nN_e : _U_nC_nN_h;
}


const double &Segment::getU_cN_cC(int state) {

    return (state == -1) ? _U_cN_cC_e : _U_cN_cC_h;
}


double Segment::getSiteEnergy(int state) {
    
  
    double result;
    if (state == -1) {
        result=getEMpoles(state) + _U_cC_nN_e;
    }
    else if (state == +1) {
        result=getEMpoles(state) + _U_cC_nN_h;
    }
     else if (state == +2) {
        result=getEMpoles(state) + _U_xX_nN_s;
    }
     else if (state == +3) {
        result=getEMpoles(state) + _U_xX_nN_t;
    }
    else {
        throw std::runtime_error(" ERROR CODE whe__00s11o__"); // blabla what do I do here?
    }
    return result;
}


void Segment::setEMpoles(int state, double energy) {

    _hasChrgState.resize(5);
    _hasChrgState[state+1] = true;
    _eMpoles[state+1] = energy;
}


double Segment::getEMpoles(int state) {

    return _eMpoles[state+1] - _eMpoles[1];
}


void Segment::AddFragment( Fragment* fragment ) {
    _fragments.push_back( fragment );
    fragment->setSegment(this);
}

void Segment::AddAtom( Atom* atom ) {
    _atoms.push_back( atom );
    atom->setSegment(this);
}

void Segment::AddPolarSite( PolarSite *pole ) {

    _polarSites.push_back(pole);
    pole->setSegment(this);

}

void Segment::AddAPolarSite( APolarSite *pole ) {

    _apolarSites.push_back(pole);
    pole->setSegment(this);

}

void Segment::calcPos() {
    vec pos = vec(0,0,0);
    double totWeight = 0.0;

    for (unsigned int i = 0; i< _atoms.size(); i++) {
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
  else if ( tag1 == "Atoms") {
    vector < Atom* > :: iterator atm;
    for (atm = _atoms.begin(); atm < _atoms.end(); ++atm) {
         int id = (*atm)->getId() % 100000;
         string name =  (*atm)->getName();
         name.resize(3);
         string resname = (*atm)->getResname();
         resname.resize(3);
         int resnr = (*atm)->getResnr();
         vec position;
         if (tag2 == "MD")      { position = (*atm)->getPos(); }
         else if (tag2 == "QM") { position = (*atm)->getQMPos(); }
         if (tag2 == "QM" && (*atm)->HasQMPart() == false) { continue; }

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
                 " ",                   // Element symbol               %2s
                 " "                    // Charge on the atom.          %2s
                 );
    }
  }
  else if ( tag1 == "Multipoles" && tag2 != "Charges") {
    vector < PolarSite* > :: iterator pol;
    for (pol = _polarSites.begin(); pol < _polarSites.end(); ++pol) {
         int id = (*pol)->getId() % 100000;
         string name =  (*pol)->getName();
         name.resize(3);
         string resname = (*pol)->getFragment()->getName();
         resname.resize(3);
         int resnr = (*pol)->getFragment()->getId() % 10000;
         vec position = (*pol)->getPos();
         
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
    else if ( tag1 == "Multipoles" && tag2 == "Charges") {
    vector < APolarSite* > :: iterator pol;
    for (pol = _apolarSites.begin(); pol < _apolarSites.end(); ++pol) {
         int id = (*pol)->getId() % 100000;
         string name =  (*pol)->getName();
         name.resize(3);
         string resname = (*pol)->getFragment()->getName();
         resname.resize(3);
         int resnr = (*pol)->getFragment()->getId() % 10000;
         vec position = (*pol)->getPos();

         fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
                      "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%4.7f\n",
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
                 " ",                   // Charge on the atom.          %2s
                 (*pol)->getQ00()
                 );
    }
  }
}


void Segment::WriteXYZ(FILE *out, bool useQMPos) {

    vector< Atom* > ::iterator ait;

    int qmatoms = 0;

    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
        if ((*ait)->HasQMPart() == true) { ++qmatoms; }
    }

    fprintf(out, "%6d \n", qmatoms);
    fprintf(out, "\n");

    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

        if ((*ait)->HasQMPart() == false && useQMPos) { continue; }

        vec     pos = (*ait)->getQMPos();
        if (!useQMPos) pos = (*ait)->getPos();
        string  name = (*ait)->getElement();

        fprintf(out, "%2s %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10);
    }


}


}}
