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


#include <votca/xtp/topology.h>
#include <votca/tools/globals.h>
#include <boost/lexical_cast.hpp>

namespace votca { namespace xtp {

Topology::Topology() : _db_id(-1), _hasPb(0), 
                       _bc(NULL), _nblist(this),
                       _isRigid(false), _isEStatified(false)  { }

// +++++++++++++++++++++ //
// Clean-Up, Destruct    //
// +++++++++++++++++++++ //

void Topology::CleanUp() {

    vector < Molecule* > ::iterator mit;
    for (mit = _molecules.begin(); mit < _molecules.end(); mit++) delete *mit;
    _molecules.clear();
    _segments.clear();
    _fragments.clear();
    _atoms.clear();
    _polarSites.clear();

    vector < SegmentType* > ::iterator sit;
    for (sit = _segmentTypes.begin(); sit < _segmentTypes.end(); sit++) {
        delete *sit;
    }
    _segmentTypes.clear();

    if (_bc) { delete(_bc); _bc = NULL; }
    _bc = new CSG::OpenBox;
    
    _nblist.Cleanup();
    _isRigid = false;
}


Topology::~Topology() {

    // clean up the list of molecules; this also deletes atoms
    vector < Molecule* > :: iterator molecule;
    for (molecule = _molecules.begin();
         molecule < _molecules.end();
         ++molecule) {
         delete *molecule;
    }
    _molecules.clear();
    _segments.clear();
    _fragments.clear();
    _atoms.clear();
    _polarSites.clear();

    vector < SegmentType* > ::iterator typeit;
    for (typeit = _segmentTypes.begin();
         typeit < _segmentTypes.end();
         typeit++) {
        delete *typeit;
    }
    _segmentTypes.clear();


}


// ++++++++ //
// Populate //
// ++++++++ //

Fragment *Topology::AddFragment(string fragment_name) {
    int fragment_id = _fragments.size() + 1;
    Fragment* fragment = new Fragment(fragment_id, fragment_name);
    _fragments.push_back(fragment);
    fragment->setTopology(this);
    return fragment;
}

Segment *Topology::AddSegment(string segment_name) {
    int segment_id = _segments.size() + 1;
    Segment* segment = new Segment(segment_id, segment_name);
    _segments.push_back(segment);
    segment->setTopology(this);
    return segment;
}

Atom *Topology::AddAtom(string atom_name) {
    int atom_id = _atoms.size() + 1;
    Atom *atom = new Atom(atom_id, atom_name);
    _atoms.push_back(atom);
    atom->setTopology(this);
    return atom;
}

PolarSite *Topology::AddPolarSite(string siteName) {
    int poleId = _polarSites.size() + 1;
    PolarSite *pole = new PolarSite(poleId, siteName);
    _polarSites.push_back(pole);
    pole->setTopology(this);
    return pole;
}

APolarSite *Topology::AddAPolarSite(string siteName) {
    int poleId = _apolarSites.size() + 1;
    APolarSite *pole = new APolarSite(poleId, siteName);
    _apolarSites.push_back(pole);
    pole->setTopology(this);
    return pole;
}

Molecule *Topology::AddMolecule(string molecule_name) {
    int molecule_id = _molecules.size() + 1;
    Molecule *molecule = new Molecule(molecule_id, molecule_name);
    _molecules.push_back(molecule);
    molecule->setTopology(this);
    return molecule;
}

SegmentType *Topology::AddSegmentType(string typeName) {
    int typeId = _segmentTypes.size() + 1;
    SegmentType *segType = new SegmentType(typeId, typeName);
    _segmentTypes.push_back(segType);
    segType->setTopology(this);
    return segType;
}



// +++++++++++++++++ //
// Periodic Boundary //
// +++++++++++++++++ //

void Topology::setBox(const matrix &box,
                      CSG::BoundaryCondition::eBoxtype boxtype) {

    // Determine box type automatically in case boxtype == typeAuto
    if(boxtype == CSG::BoundaryCondition::typeAuto) {
        boxtype = AutoDetectBoxType(box);
    }

    if(_hasPb) {
        if (votca::tools::globals::verbose) {
            cout << "Removing periodic box. Creating new... " << endl;
        }
        delete _bc;
    }

    switch(boxtype) {
       case CSG::BoundaryCondition::typeTriclinic:
            _bc = new CSG::TriclinicBox();
            break;
       case CSG::BoundaryCondition::typeOrthorhombic:
            _bc = new CSG::OrthorhombicBox();
            break;
       default:
            _bc = new CSG::OpenBox();
            break;
       }

    _bc->setBox(box);
    _hasPb = true;
}


CSG::BoundaryCondition::eBoxtype Topology::AutoDetectBoxType(const matrix &box){

    // Set box type to OpenBox in case "box" is the zero matrix,
    // to OrthorhombicBox in case "box" is a diagonal matrix,
    // or to TriclinicBox otherwise

    if(box.get(0,0)==0 && box.get(0,1)==0 && box.get(0,2)==0 &&
       box.get(1,0)==0 && box.get(1,1)==0 && box.get(1,2)==0 &&
       box.get(2,0)==0 && box.get(2,1)==0 && box.get(2,2)==0) {

            cout << "WARNING: No box vectors specified in trajectory."
               "Using open-box boundary conditions. " << endl;
            return CSG::BoundaryCondition::typeOpen;
    }

    else if(box.get(0,1)==0 && box.get(0,2)==0 &&
            box.get(1,0)==0 && box.get(1,2)==0 &&
            box.get(2,0)==0 && box.get(2,1)==0) {

            return CSG::BoundaryCondition::typeOrthorhombic;
    }

    else {
        return CSG::BoundaryCondition::typeTriclinic;
    }

    return CSG::BoundaryCondition::typeOpen;
}


vec Topology::PbShortestConnect(const vec &r1, const vec &r2) const {
    return _bc->BCShortestConnection(r1, r2);
}



bool Topology::Rigidify() {

    if (!_canRigidify) {
        cout << endl
             << "... ... ABORT: Request to rigidify system, but no QM "
                "coordinates provided. ";
        return 0;
    }
    else {

        // Rigidify segments
        vector<Segment*> ::iterator sit;
        for (sit = _segments.begin();
             sit < _segments.end();
             sit++) {

             (*sit)->Rigidify();
        }

        cout << endl
             << "... ... Rigidified " << _segments.size() << " segments. "
             << flush;

        if (this->NBList().size() > 0) {

        // Rigidify pairs
        // [ Why this is needed: Orientation matrices for fragments are not
        //   not written to the state file, and hence lost whenever a frame
        //   is saved and reloaded. When reading in a new frame from the
        //   database, pairs are created from scratch based on two segment
        //   IDs each. Then the pair constructor is called to check whether
        //   the pair is formed across the periodic boundary. If so, it
        //   creates a ghost from the second partner. This ghost is a new
        //   segment which is just accessible from within the pair, but not
        //   from the topology; i.e. it is not stored in any segment containers.
        //   Since at this point, fragments do not store rotation matrices yet,
        //   the ghost - very much deconnected from its originator - does not,
        //   either. Therefore it should not be forgotten here. --- A way out
        //   would be to rigidify the topology within StateSaver::ReadFrame,
        //   after atoms have been created, but before pairs are created. ]

            QMNBList &nblist = this->NBList();

            QMNBList::iterator pit;
            int count = 0;
            for (pit = nblist.begin(); pit != nblist.end(); pit++) {

                QMPair *qmpair = *pit;
                if (qmpair->HasGhost()) {
                    count++;



                    qmpair->Seg2PbCopy()->Rigidify();
                }
            }

          cout << endl
               << "... ... Rigidified " << count << " ghosts. "
               << flush;
        }

        _isRigid = true;
        return 1;
    }
}







void Topology::PrintInfo(ostream &out) {
        cout << endl;

        cout << "MD|QM Topology info "
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        if (_bc) {
        cout << "Periodic Box:        "  << this->getBox().get(0,0)
                                << " "   << this->getBox().get(0,1)
                                << " "   << this->getBox().get(0,2)
                                << " | " << this->getBox().get(1,0)
                                << " "   << this->getBox().get(1,1)
                                << " "   << this->getBox().get(1,2)
                                << " | " << this->getBox().get(2,0)
                                << " "   << this->getBox().get(2,1)
                                << " "   << this->getBox().get(2,2)
                                << endl;
        }
        else {
        cout << "No periodic boundary specified. " << endl;
        }
        cout << "Database ID:         " << this->getDatabaseId() << endl;
        cout << "Step number:         " << this->getStep() << endl;
        cout << "Time:                " << this->getTime() << endl;
        cout << "# Molecules          " << this->Molecules().size() << endl;
        cout << "# Segments           " << this->Segments().size() << endl;
        cout << "# Fragments          " << this->Fragments().size() << endl;
        cout << "# Atoms              " << this->Atoms().size() << endl;
        cout << "# Pairs              " << this->NBList().size() << endl;
        cout << endl;
}



void Topology::PrintInfo(FILE *out) {

    fprintf(out, "Topology Database ID %3d \n", this->getDatabaseId());
    fprintf(out, "Periodic Box: %2.4f %2.4f %2.4f | %2.4f %2.4f %2.4f | %2.4f %2.4f %2.4f \n",
            this->getBox().get(0,0),
            this->getBox().get(0,1),
            this->getBox().get(0,2),
            this->getBox().get(1,0),
            this->getBox().get(1,1),
            this->getBox().get(1,2),
            this->getBox().get(2,0),
            this->getBox().get(2,1),
            this->getBox().get(2,2) );

    fprintf(out, "\tStep number %7d \n", this->getStep());
    fprintf(out, "\tTime        %2.4f \n", this->getTime());
    fprintf(out, "\t# Molecules %7ld \n", this->Molecules().size());
    fprintf(out, "\t# Segments  %7ld \n", this->Segments().size());
    fprintf(out, "\t# Atoms     %7ld \n", this->Atoms().size());
    fprintf(out, "\t# Pairs     %7ld \n", this->NBList().size());

}


void Topology::WritePDB(FILE *out, string tag) {

    if (tag == "segments") {

    vector < Segment* > :: iterator seg;
    for (seg = _segments.begin(); seg < _segments.end(); ++seg){
         int id = (*seg)->getId();
         string name =  (*seg)->getName();
         name.resize(3);
         string resname = (*seg)->getMolecule()->getName();
         resname.resize(3);
         int resnr = (*seg)->getMolecule()->getId();
         vec position = (*seg)->getPos();

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
    } /* exit tag segments */

    else {
        ;
    }
}



}}
