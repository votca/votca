/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include "Md2QmEngine.h"
#include <votca/csg/boundarycondition.h>
#include <votca/tools/globals.h>

/**
 * Clears all engine template ('type') containers.
 */
Md2QmEngine::~Md2QmEngine() {

    vector < CTP::Molecule* > :: iterator molecule;

    // Clean up list of molecule types
    for (molecule = _molecule_types.begin();
         molecule < _molecule_types.end(); ++molecule){
         delete *molecule;
    }
    _molecule_types.clear();

    // clean up maps
    _map_MoleculeName_MoleculeType.clear();
    _map_MoleculeMDName_MoleculeName.clear();
}



/**
 * Reads in mapping .xml file to generate molecule, segment, fragment and
 * atom templates ('types') that contain all necessary information to
 * start mapping, except coordinates.
 * @param xmlfile
 */
void Md2QmEngine::Initialize(const string &xmlfile) {

    Property typology;
    load_property_from_xml(typology, xmlfile.c_str());
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // XML to Types: Molecules => Segments => Fragments => Atoms //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    string key = "topology.molecules.molecule";
    list<Property *> molecules = typology.Select(key);
    list<Property *>::iterator it_molecule;
    int molecule_id = 1;
    int qmunit_id  = 1; // Counts segment types

    for ( it_molecule = molecules.begin();
          it_molecule != molecules.end();
          ++it_molecule ) {

       CTP::Molecule *molecule = AddMoleculeType(molecule_id++, *it_molecule);
       string molMdName = (*it_molecule)->get("mdname").as<string>();

       // +++++++++++++ //
       // Load segments //
       // +++++++++++++ //

       key = "segments.segment";
       list<Property *> segments = (*it_molecule)->Select(key);
       list<Property *>::iterator it_segment;
       int segment_id = 1;       
       int md_atom_id = 1; // <- atom id count with respect to molecule

       for ( it_segment = segments.begin();
             it_segment != segments.end();
             ++it_segment ) {

         // Create new segment + associated type (QM Unit)
         CTP::Segment *segment = AddSegmentType(segment_id++, *it_segment );
         CTP::SegmentType *qmUnit = AddQMUnit(qmunit_id, *it_segment );
         ++qmunit_id;

         segment->setType(qmUnit);
         molecule->AddSegment(segment);

         // Load internal (i.e. QM-) coord.s and MOO-related properties
         string qmcoordsFile = "";
         map<int, pair<string, vec> > intCoords;
         if ( (*it_segment)->exists("qmcoords") ) {
            qmcoordsFile = (*it_segment)->get("qmcoords").as<string>();            
            //  QM ID    Element   Position
            this->getIntCoords(qmcoordsFile, intCoords);
         }
         
         // ++++++++++++++ //
         // Load fragments //
         // ++++++++++++++ //
         
         map<string,bool> fragname_isTaken;
         
         key = "fragments.fragment";
         list<Property *> fragments = (*it_segment)->Select(key);
         list<Property *>::iterator it_fragment;
         int fragment_id = 1;

         for ( it_fragment = fragments.begin();
               it_fragment != fragments.end();
               ++it_fragment ) {

            // Create new fragment
            CTP::Fragment* fragment=AddFragmentType(fragment_id++,*it_fragment);
            segment->AddFragment( fragment );
            
            // Verify that this fragment name is not taken already
            try {
                // This should throw
                bool taken = fragname_isTaken.at(fragment->getName());
                cout << "ERROR Fragment name '" << fragment->getName()
                     << "' in segment '" << segment->getName()
                     << "' occurs more than once." << endl;
                if (taken)
                    throw runtime_error("(see above, naming collision)");
            }
            catch (const std::exception& out_of_range) {
                fragname_isTaken[fragment->getName()] = true;
            }
            
            // Load local-frame definition
            vector<int> trihedron = (*it_fragment)->get("localframe")
                                                .as< vector<int> >();
            while (trihedron.size() < 3) {
                trihedron.push_back(-1);
            }
            fragment->setTrihedron(trihedron);


             // ++++++++++ //
             // Load atoms //
             // ++++++++++ //

             string mdatoms = (*it_fragment)->get("mdatoms").as<string>();
             string qmatoms = (*it_fragment)->get("qmatoms").as<string>();
             string weights = (*it_fragment)->get("weights").as<string>();

             Tokenizer tok_md_atoms(mdatoms, " \t\n");
             Tokenizer tok_qm_atoms(qmatoms, " \t\n");
             Tokenizer tok_weights(weights, " \t\n");

             vector <string> md_atoms_info;
             vector <string> qm_atoms_info;
             vector <string> atom_weights;;

             tok_md_atoms.ToVector(md_atoms_info);
             tok_qm_atoms.ToVector(qm_atoms_info);
             tok_weights.ToVector(atom_weights);

             //assert(md_atoms_info.size() == qm_atoms_info.size());
             //assert(md_atoms_info.size() == atom_weights.size());

             if ( (md_atoms_info.size() != qm_atoms_info.size()) ||
                  (md_atoms_info.size() != atom_weights.size() ) ) {
                 cout << "ERROR: "
                      << "Could not allocate MD atoms to QM atoms or weights"
                      << " in fragment " << fragment->getName()
                      << " in segment " << segment->getName()
                      << " in molecule " << molMdName
                      << " due to inconsistent number of columns"
                      << " (MD: " << md_atoms_info.size() << ","
                      << " QM: " << qm_atoms_info.size() << ")"
                      << " Weights: " << atom_weights.size() << ")."
                      << endl;
                 cout << "NOTE: "
                      << "To define an MD atom without QM counterpart, insert "
                      << "a single ':' in the associated QM-atoms column and "
                      << "specify a mapping weight of 0."
                      << endl;
                 throw runtime_error( "Inconsistency in mapping file." );
             }

             vector<string> ::iterator it_md_atom_name;
             vector<string> ::iterator it_qm_atom_name;
             vector<string> ::iterator it_atom_weight;

             for ( it_md_atom_name = md_atoms_info.begin(),
                   it_qm_atom_name = qm_atoms_info.begin(),
                   it_atom_weight  = atom_weights.begin();
                   it_md_atom_name != md_atoms_info.end();
                   ++it_md_atom_name,
                   ++it_qm_atom_name,
                   ++ it_atom_weight) {

               // ++++++++++++++++++ //
               // Create single atom //
               // ++++++++++++++++++ //

               Tokenizer tok_md((*it_md_atom_name), ":");
               Tokenizer tok_qm((*it_qm_atom_name), ":");

               vector<string> md_atom_specs;
               vector<string> qm_atom_specs;

               tok_md.ToVector( md_atom_specs );
               tok_qm.ToVector( qm_atom_specs );

               // MD atom information
               int residue_number  = boost::lexical_cast<int>(md_atom_specs[0]);
               string residue_name = md_atom_specs[1];
               string md_atom_name = md_atom_specs[2];
               
               // QM atom information
               bool hasQMPart = false;
               int qm_atom_id = -1;
               vec qmPos = vec(0,0,0);
               // The atomic element is first taken as first character of
               // the MD name; for atoms with QM counterpart, the element
               // is read from the coordinates file
               string element = md_atom_name.substr(0,1);
               // Check whether MD atom has QM counterpart
               if (qm_atom_specs.size() == 2) {
                   hasQMPart = true;
                   qm_atom_id = boost::lexical_cast<int>(qm_atom_specs[0]);

                   // Look up qm coordinates in table created from xyz file
                   try {
                       qmPos = intCoords.at(qm_atom_id).second;
                       element = intCoords.at(qm_atom_id).first;
                       // Check whether elements of md and qm match
                       if (intCoords.at(qm_atom_id).first.substr(0,1)
                           != md_atom_name.substr(0,1) ) {
                           cout << "WARNING: Atom " <<md_atom_name << " in mol. "
                                << molMdName << " appears to have element type "
                                << md_atom_name.substr(0,1)
                                << ", but QM partner (ID " << qm_atom_id
                                << ") has element type "
                                << intCoords.at(qm_atom_id).first << endl;
                       }
                   }
                   catch (const std::exception& out_of_range) {
                       ; // No QM coordinates specified
                   }
               }
               
               // Mapping weight
               double weight = boost::lexical_cast<double>(*it_atom_weight);
               if (!hasQMPart && weight != 0) {
                   cout << "ERROR: "
                        << "Atom " << md_atom_name << " in residue "
                        << residue_name << " in molecule " << molMdName
                        << " has no QM counterpart despite non-zero weight. "
                        << endl;
                   throw runtime_error( "Error in mapping file" );
               }

               // Create atom
               CTP::Atom *atom = AddAtomType(molecule,
                                             residue_name, residue_number,
                                             md_atom_name, md_atom_id++,
                                             hasQMPart,    qm_atom_id,
                                             qmPos,        element,
                                             weight);

               try {
                   this->_map_mol_resNr_atm_atmType.at(molMdName)
                                                   .at(residue_number)
                                                   .at(md_atom_name);
                   // If this succeeded, atom has already been defined:
                   cout << "ERROR: "
                        << "Ambiguous atom definition in molecule "
                        << molMdName << ": "
                        << "Atom " << md_atom_name << " in residue "
                        << residue_number << " exists more than once. "
                        << endl;
                   throw runtime_error( "Ambiguity in atom definition." );
               }
               catch ( const std::exception& out_of_range ) {
                   this->_map_mol_resNr_atm_atmType[molMdName]
                                                   [residue_number]
                                                   [md_atom_name] = atom;
               }

               fragment->AddAtom(atom);
               segment->AddAtom(atom);
               molecule->AddAtom(atom);

             } /* exit loop over atoms */

          } /* exit loop over fragments */

       } /* exit loop over segments */

    } /* exit loop over molecules */
}



/**
 * Using the templates for molecules, segments, ... generated in
 * ::Initialize, Md2Qm maps the MD topology to the initially empty
 * MD/QM Topology. Calls ::MoleculeFactory and ::ExportMolecule.
 * @param mdtop
 * @param qmtop
 */
void Md2QmEngine::Md2Qm(CSG::Topology *mdtop, CTP::Topology *qmtop) {

    qmtop->CleanUp();

    // Create periodic box
    qmtop->setBox(mdtop->getBox());

    // Set trajectory meta data
    qmtop->setStep(mdtop->getStep());
    qmtop->setTime(mdtop->getTime());
    qmtop->setCanRigidify(true);

    // Add types (=> Segment types / QM units)
    vector< CTP::SegmentType* > ::iterator typeit;
    for (typeit = this->_qmUnits.begin();
         typeit != this->_qmUnits.end();
         typeit++) {

        string name = (*typeit)->getName();
        string basis = (*typeit)->getBasisName();
        string orbitals = (*typeit)->getOrbitalsFile();
        string qmcoords = (*typeit)->getQMCoordsFile();
        bool canRigidify = (*typeit)->canRigidify();
        CTP::SegmentType *segType = qmtop->AddSegmentType(name);
        segType->setBasisName(basis);
        segType->setOrbitalsFile(orbitals);
        segType->setQMCoordsFile(qmcoords);
        segType->setCanRigidify(canRigidify);
        if (!canRigidify) { qmtop->setCanRigidify(false); }
    }

    // Populate topology in a trickle-down manner
    // (i.e. molecules => ... ... => atoms)
    CSG::MoleculeContainer::iterator mit;
    for (mit = mdtop->Molecules().begin();
         mit != mdtop->Molecules().end();
         mit++ ) {

         // MD molecule + name
         CSG::Molecule *molMD = *mit;
         string nameMolMD = molMD->getName();

         // Find QM counterpart
         CTP::Molecule *molQM = this->MoleculeFactory(molMD);
         string nameMolQM = molQM->getName();

         // Generate and export
         //CTP::Molecule *product = 
	 (void)this->ExportMolecule(molQM, qmtop);
    }
}



/**
 * Takes MD molecule, finds QM counterpart, fills out atom positions.
 * @param molMDTemplate
 * @return molQMTemplate
 */
CTP::Molecule *Md2QmEngine::MoleculeFactory(CSG::Molecule *molMDTemplate) {

    string nameMolQM = this->getMoleculeName(molMDTemplate->getName());
    CTP::Molecule *molQMTemplate = this->getMoleculeType(nameMolQM);

    int resnrOffset = molMDTemplate->getBead(0)->getResnr();

    for (int i = 0; i < molMDTemplate->BeadCount(); i++) {
        CSG::Bead *atomMD = molMDTemplate->getBead(i);
        try {
            CTP::Atom *counterpart =
                 this->getAtomType(molMDTemplate->getName(),
                                   atomMD->getResnr()-resnrOffset+1,
                                   atomMD->getName());
            counterpart->setPos(atomMD->getPos());
        }
        catch (const std::exception& out_of_range) {
            cout << "WARNING: No mapping instruction found for atom "
                 << atomMD->getName() << " in residue number "
                 << atomMD->getResnr() << " in molecule "
                 << molMDTemplate->getName() << ". Skipping..."
                 << endl;
        }
    }
    return molQMTemplate;
}



/**
 * Takes QM molecule template as reference, and step by step adds 'new'
 * copies of all segments, fragments, atoms within that molecule to
 * the target topology.
 * @param refMol
 * @param qmtop
 * @return
 */
CTP::Molecule *Md2QmEngine::ExportMolecule(CTP::Molecule *refMol,
                                           CTP::Topology *qmtop) {

    CTP::Molecule *newMol = qmtop->AddMolecule(refMol->getName());

    // Note: The topology is responsible for allocating IDs to
    //       molecules, segments, ...
    
    vector<CTP::Segment *> ::iterator segIt;
    for (segIt = refMol->Segments().begin();
         segIt < refMol->Segments().end();
         segIt++) {

        CTP::Segment *refSeg = *segIt;
        CTP::Segment *newSeg = qmtop->AddSegment(refSeg->getName());
        newSeg->setType( qmtop->getSegmentType(refSeg->getType()->getId()) );

        vector<CTP::Fragment *> ::iterator fragIt;
        for (fragIt = refSeg->Fragments().begin();
             fragIt < refSeg->Fragments().end();
             fragIt++) {

            CTP::Fragment *refFrag = *fragIt;
            CTP::Fragment *newFrag = qmtop->AddFragment(refFrag->getName());
            newFrag->setTrihedron(refFrag->getTrihedron());

            vector<CTP::Atom *> ::iterator atomIt;
            for (atomIt = refFrag->Atoms().begin();
                 atomIt < refFrag->Atoms().end();
                 atomIt++) {

                CTP::Atom *refAtom = *atomIt;
                CTP::Atom *newAtom = qmtop->AddAtom(refAtom->getName());

                newAtom->setWeight(refAtom->getWeight());
                newAtom->setResnr(refAtom->getResnr());
                newAtom->setResname(refAtom->getResname());
                newAtom->setElement(refAtom->getElement());
                if (refAtom->HasQMPart()) {
                    newAtom->setQMPart(refAtom->getQMId(),
                                       refAtom->getQMPos());
                }
                newAtom->setPos(refAtom->getPos());

                newFrag->AddAtom(newAtom);
                newSeg->AddAtom(newAtom);
                newMol->AddAtom(newAtom);
            } /* exit loop over template atoms */
            newFrag->calcPos();
            newSeg->AddFragment(newFrag);
            newMol->AddFragment(newFrag);
        } /* exit loop over template fragments */
        newSeg->calcPos();
        newMol->AddSegment(newSeg);
    } /* exit loop over template molecules */

    return newMol;
}



CTP::Atom *Md2QmEngine::AddAtomType(CTP::Molecule *owner,
                                    string residue_name,  int residue_number,
                                    string md_atom_name,  int md_atom_id,
                                    bool hasQMPart,       int qm_atom_id,
                                    vec qmPos,            string element,
                                    double weight) {
    CTP::Atom* atom = new CTP::Atom(owner,
                                    residue_name,         residue_number,
                                    md_atom_name,         md_atom_id,
                                    hasQMPart,            qm_atom_id,
                                    qmPos,                element,
                                    weight);
    _atom_types.push_back(atom);
    return atom;
}

CTP::Fragment *Md2QmEngine::AddFragmentType(int fragment_id,
                                            Property *property) {
    string fragment_name = property->get("name").as<string>();
    CTP::Fragment* fragment = new CTP::Fragment(fragment_id, fragment_name);
    _fragment_types.push_back(fragment);
    return fragment;
}

CTP::Segment  *Md2QmEngine::AddSegmentType(int segment_id,
                                           Property *property) {
    string segment_name = property->get("name").as<string>();
    CTP::Segment* segment = new CTP::Segment(segment_id, segment_name);
    _segment_types.push_back(segment);
    return segment;
}

CTP::Molecule *Md2QmEngine::AddMoleculeType(int molecule_id,
                                            Property *property) {
    string molecule_name = property->get("name").as<string>();
    string molecule_mdname = property->get("mdname").as<string>();

    CTP::Molecule *molecule = new CTP::Molecule(molecule_id, molecule_name);
    _molecule_types.push_back(molecule);
    // TODO check if the name is already there
    _map_MoleculeMDName_MoleculeName[molecule_mdname] = molecule_name;
    _map_MoleculeName_MoleculeType[molecule_name] = molecule;
    return molecule;
}

CTP::SegmentType *Md2QmEngine::AddQMUnit(int qmunit_id, Property *property) {
    string qmCoordsFile = "nofile";
    string orbitalsFile = "nofile";
    string basisSetName = "noname";
    vector<int> torbNrs;

    bool canRigidify = false;
        
    string qmunit_name = property->get("name").as<string>();
    
    if (property->exists("qmcoords")) {
        qmCoordsFile = property->get("qmcoords").as<string>();
        canRigidify = true;
    }
    if (property->exists("orbitals")) {
        orbitalsFile = property->get("orbitals").as<string>();
    }
    if (property->exists("basisset")) {
        basisSetName = property->get("basisset").as<string>();
    }

    CTP::SegmentType* qmUnit = new CTP::SegmentType(qmunit_id,    qmunit_name,
                                                    basisSetName, orbitalsFile,
                                                    qmCoordsFile, canRigidify);
    _qmUnits.push_back(qmUnit);
    return qmUnit;
}


CTP::Molecule *Md2QmEngine::getMoleculeType(const string &name) {
    try {
        return _map_MoleculeName_MoleculeType.at(name);
    }
    catch ( const std::exception& out_of_range ) {
        cout << "WARNING: Molecule '" << name
             << "' not included in mapping definition. Skipping... ";
        cout << endl;
        return NULL;
    }    
}

const string &Md2QmEngine::getMoleculeName(const string &mdname) {
    try {
        return _map_MoleculeMDName_MoleculeName.at(mdname);
    }
    catch ( const std::exception& out_of_range ) {
        return mdname;
    }
}

CTP::Atom *Md2QmEngine::getAtomType(const string &molMdName,
                                    int resNr, const string &mdAtomName) {
    return this->_map_mol_resNr_atm_atmType.at(molMdName)
                                           .at(resNr)
                                           .at(mdAtomName);
}

void Md2QmEngine::getIntCoords(string &file,
                               map<int, pair<string,vec> > &intCoords) {


    int atomCount = 0;

    std::string line;
    std::ifstream intt;
    intt.open(file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {
            std::getline(intt, line);

            vector< string > split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);
            if ( !split.size()      ||
                  split.size() != 4 ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

            // Interesting information written here: e.g. 'C 0.000 0.000 0.000'
            atomCount++;
            string element = split[0];
            double x = boost::lexical_cast<double>( split[1] ) / 10.; //°A to NM
            double y = boost::lexical_cast<double>( split[2] ) / 10.;
            double z = boost::lexical_cast<double>( split[3] ) / 10.;
            vec qmPos = vec(x,y,z);

            pair<string, vec> qmTypePos(element, qmPos);
            intCoords[atomCount] = qmTypePos;
        }
    }
    else {
        throw std::runtime_error("No such file: '"+file+"'.");
    }
}


void Md2QmEngine::PrintInfo() {

    vector<CTP::Molecule*>::iterator mit;
    vector<CTP::Segment*>::iterator sit;
    vector<CTP::Fragment*>::iterator fit;

    cout << "Summary ~~~~~"
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Created "
         << _molecule_types.size() << " molecule type(s): ";
    for (mit = _molecule_types.begin();
         mit != _molecule_types.end();
         mit++) {
         cout << "[ " << (*mit)->getName() << " ]";
    }
    cout << endl <<  "with a total of "
         << _segment_types.size()  << " segment(s), "
         << _fragment_types.size() << " fragments, "
         << _atom_types.size() << " atoms. \n" << endl;

    map < string, string > ::iterator mssit;
    for (mssit = this->_map_MoleculeMDName_MoleculeName.begin();
         mssit != this->_map_MoleculeMDName_MoleculeName.end();
         mssit ++) {
         cout << "MD [ " << mssit->first << " ] mapped to "
              << "QM [ " << mssit->second << " ] \n" << endl;
    }

    for (mit = _molecule_types.begin();
         mit != _molecule_types.end();
         mit++) {
         CTP::Molecule *mol = *mit;
         cout << "[ " << mol->getName() << " ]" << endl;

         for (sit = mol->Segments().begin();
              sit != mol->Segments().end();
              sit++) {
              CTP::Segment *seg = *sit;
              cout << " - Segment [ " << seg->getName() << " ]"
                   << " ID " << seg->getId() << endl;

              for (fit = seg->Fragments().begin();
                   fit != seg->Fragments().end();
                   fit++ ) {
                   CTP::Fragment *frag = *fit;
                   cout << "   - Fragment [ " << frag->getName() << " ]"
                        << " ID " << frag->getId() << ": "
                        << frag->Atoms().size() << " atoms " << endl;
              }
         }
    }

    if (! votca::tools::globals::verbose ) { return; }

    cout << endl << "Mapping table"
                    " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                 << endl;
    map < string, map < int, map < string, CTP::Atom* > > > ::iterator it0;
    for (it0 = this->_map_mol_resNr_atm_atmType.begin();
         it0 != this->_map_mol_resNr_atm_atmType.end();
         it0 ++) {
         map < int, map < string, CTP::Atom* > > ::iterator it1;
         for (it1 = it0->second.begin();
              it1 != it0->second.end();
              it1++) {
              map < string, CTP::Atom* > ::iterator it2;
              for (it2 = it1->second.begin();
                   it2 != it1->second.end();
                   it2++) {

       printf( "MD Molecule %4s | Residue %2d | Atom %3s "
                               "| ID %3d => QM ID %3d \n",
               it0->first.c_str(),
               it2->second->getResnr(),
               it2->first.c_str(),
               it2->second->getId(),
               it2->second->getQMId());
             }
         }
    }
}

void Md2QmEngine::CheckProduct(CTP::Topology *outtop, const string &pdbfile) {

    string md_pdb = "md_" + pdbfile;
    string qm1_pdb = "qm_rigid_" + pdbfile;
    string qm2_pdb = "qm_conjg_" + pdbfile;
    FILE *outPDB = NULL;

    // Atomistic PDB
    outPDB = fopen(md_pdb.c_str(), "w");

    vector<CTP::Molecule*> ::iterator molIt;
    for (molIt = outtop->Molecules().begin();
         molIt < outtop->Molecules().end();
         molIt++) {
        CTP::Molecule *mol = *molIt;
        mol->WritePDB(outPDB);
    }

    fprintf(outPDB, "\n");
    fclose(outPDB);


    // Fragment PDB
    outPDB = fopen(qm1_pdb.c_str(), "w");

    vector<CTP::Segment*> ::iterator segIt;
    for (segIt = outtop->Segments().begin();
         segIt < outtop->Segments().end();
         segIt++) {
        CTP::Segment *seg = *segIt;
        seg->WritePDB(outPDB);
    }

    fprintf(outPDB, "\n");
    fclose(outPDB);


    // Segment PDB
    outPDB = fopen(qm2_pdb.c_str(), "w");
    outtop->WritePDB(outPDB);
    fprintf(outPDB, "\n");
    fclose(outPDB);


    if (TOOLS::globals::verbose) {
        cout << endl;
        this->PrintInfo();
        cout << endl;
        outtop->PrintInfo(cout);
    }
}

