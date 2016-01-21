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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <votca/xtp/molecule.h>
#include <votca/xtp/atom.h>

using namespace std;

namespace votca { namespace xtp {

/// Destructor
Molecule::~Molecule() {

    vector < Segment* > ::iterator segit;
    for (segit = this->Segments().begin();
            segit < this->Segments().end();
            segit ++) {
        delete *segit;
    }
    _segments.clear();
    _fragments.clear();
    _atoms.clear();
    
}

/// Returns ID of the molecule
const int &Molecule::getId() {
        return _id;
}

/// Returns the name of the molecule
const string &Molecule::getName() {
        return _name;
}

void Molecule::AddSegment( Segment* segment ) {
    _segments.push_back( segment );
    segment->setMolecule(this);
}

void Molecule::AddFragment(Fragment* fragment ) {
    _fragments.push_back( fragment );
    fragment->setMolecule(this);
}

void Molecule::AddAtom(Atom *atom ) {
    _atoms.push_back(atom);
    atom->setMolecule(this);
}

/// Returns a pointer to the atom
Atom *Molecule::getAtom(const int &id) {
    throw runtime_error( string("Not implemented") ); 
}

/// Returns atom type
const string &Molecule::getAtomType(const int &id) {
    throw runtime_error( string("Not implemented") ); 
}

/// Returns atom position
const vec Molecule::getAtomPosition(const int &id) {
    throw runtime_error( string("Not implemented") ); 
}

/// Returns number of atoms in the molecule
int Molecule::NumberOfAtoms() {
    return _atoms.size(); 
}

/// Read atom types and coordinates from a file
void Molecule::ReadXYZ(string filename)
{
    ifstream in;
    double x,y,z;
    //int natoms, id;
    string label, type;
    vec pos;
    
    cout << " Reading molecular coordinates from " << filename << endl;
    in.open(filename.c_str(), ios::in );
    if ( !in ) throw runtime_error( string("Error reading coordinates from: ") 
                                        + filename); 
    
    //id = 1;
    while ( in.good() ) { // keep reading until end-of-file
        in >> type;
        in >> x;
        in >> y;
        in >> z;
        if( in.eof() ) break;
        //cout << type << ":" << x << ":" << y << ":" << z << endl;
        
        // creating atoms and adding them to the molecule
        //Atom *pAtom = new Atom( this, id++, type, 1, 0.12, 0.12);
        //vec position(x,y,z);
        //pAtom->setPos( position );                       
        //AddAtom( pAtom );  

    }      
    in.close();        
}

/// Writes a PDB file
void Molecule::WritePDB( FILE *out ) {
    // out.setf(ios::fixed);
    
    vector < Atom* > :: iterator atom;
    for (atom = _atoms.begin(); atom < _atoms.end(); ++atom){
         int id = (*atom)->getId();      
         string name =  (*atom)->getName();
         string resname = (*atom)->getResname();
         int resnr = (*atom)->getResnr();
         vec position = (*atom)->getPos();  
         
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

/*
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.00           H

 1 -  6        Record name     "ATOM  "                                            
 7 - 11        Integer         Atom serial number.                   
13 - 16        Atom            Atom name.                            
17             Character       Alternate location indicator.         
18 - 20        Residue name    Residue name.                         
22             Character       Chain identifier.                     
23 - 26        Integer         Residue sequence number.              
27             AChar           Code for insertion of residues.       
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
55 - 60        Real(6.2)       Occupancy.                            
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
73 - 76        LString(4)      Segment identifier, left-justified.   
77 - 78        LString(2)      Element symbol, right-justified.      
79 - 80        LString(2)      Charge on the atom.       
*/

}}
