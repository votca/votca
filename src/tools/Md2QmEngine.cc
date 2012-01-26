#include "Md2QmEngine.h"



Md2QmEngine::~Md2QmEngine() {

    vector < CTP::Fragment* > :: iterator fragment;
    vector < CTP::Segment* > :: iterator segment;
    vector < CTP::Molecule* > :: iterator molecule;

    // Clean up list of fragment types
    // cout << " * Typology: deleting fragment types: ";
    for (fragment = _fragment_types.begin();
         fragment < _fragment_types.end();
         ++fragment) {
         //cout << (*fragment)->getId() << ":"
         //     << (*fragment)->getName() << " ";
         delete *fragment;
    }
    //cout << endl;
    _fragment_types.clear();

    // Clean up list of segment types
    // cout << " * Typology: deleting segment types: ";
    for (segment = _segment_types.begin();
         segment < _segment_types.end();
         ++segment) {
         //cout << (*segment)->getId() << ":"
         //     << (*segment)->getName() << " ";
         delete *segment;
    }
    //cout << endl;
    _segment_types.clear();

    // Clean up list of molecule types; this also deletes atoms
    // cout << " * Typology: deleting molecule types: ";
    for (molecule = _molecule_types.begin();
         molecule < _molecule_types.end(); ++molecule){
         //cout << (*molecule)->getId() << ":"
         //     << (*molecule)->getName() << " ";
         delete *molecule;
    }
    //cout << endl;
    _molecule_types.clear();

    // clean up the map of the molecule md name : pointer to molecule type
    _map_MoleculeName_MoleculeType.clear();
    _map_MoleculeMDName_MoleculeName.clear();

}


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

       for ( it_segment = segments.begin();
             it_segment != segments.end();
             ++it_segment ) {

         // Create a new segment
         CTP::Segment *segment = AddSegmentType(segment_id++, *it_segment );
         molecule->AddSegment( segment );

         // ++++++++++++++ //
         // Load fragments //
         // ++++++++++++++ //

         key = "fragments.fragment";
         list<Property *> fragments = (*it_segment)->Select(key);
         list<Property *>::iterator it_fragment;
         int fragment_id = 1;

         for ( it_fragment = fragments.begin();
               it_fragment != fragments.end();
               ++it_fragment ) {

             // Create new fragment
             CTP::Fragment* fragment=AddFragmentType(fragment_id,*it_fragment);
             segment->AddFragment( fragment );

             // ++++++++++ //
             // Load atoms //
             // ++++++++++ //

             string mdatoms = (*it_fragment)->get("mdatoms").as<string>();
             string qmatoms = (*it_fragment)->get("qmatoms").as<string>();
             string weights = (*it_fragment)->get("weights").as<string>();

             Tokenizer tok_md_atoms(mdatoms, " ");
             Tokenizer tok_qm_atoms(qmatoms, " ");
             Tokenizer tok_weights(weights, " ");

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
                      << " in segment " << segment->getName()
                      << " in molecule " << molMdName
                      << " due to inconsistent number of columns."
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
             int md_atom_id = 1;

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
               // Check whether MD atom has QM counterpart
               if (qm_atom_specs.size() == 2) {
                   hasQMPart = true;
                   qm_atom_id = boost::lexical_cast<int>(qm_atom_specs[0]);
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
               CTP::Atom *atom = AddAtomType(molecule,     residue_number,
                                             md_atom_name, md_atom_id,
                                             hasQMPart,    qm_atom_id,
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
               catch ( out_of_range ) {
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







void Md2QmEngine::Md2Qm(CSG::Topology *mdtop, CTP::Topology *qmtop) {

    int globalMolID = 0;
    int globalSegID = 0;
    int globalFrgID = 0;
    int globalAtmID = 0;

    CSG::MoleculeContainer::iterator mit;
    for (mit = mdtop->Molecules().begin();
         mit != mdtop->Molecules().end();
         mit++ ) {

         // MD molecule + name
         CSG::Molecule *molMD = *mit;
         string nameMolMD = molMD->getName();

         // Find QM counterpart
         string nameMolQM = this->getMoleculeName(nameMolMD);
         CTP::Molecule *MolQM = this->getMoleculeType(nameMolQM);

         int resnrOffset = molMD->getBead(0)->getResnr();




        for (int i = 0; i < molMD->BeadCount(); i++) {
             CSG::Bead* atomMD = molMD->getBead(i);







/*
            cout << "ATOM getId  " << atomMD->getId();
            cout << " | getName  " << atomMD->getName();
            cout << " | getMol   " << atomMD->getMolecule()->getId();
            cout << " | MolName  " << atomMD->getMolecule()->getName();
            cout << " | getResnr " << atomMD->getResnr();
            cout << " | getType  " << atomMD->getType()->getName();
            cout << " | getPos   " << atomMD->getPos();
            cout << endl;*/
        }
    }
}















CTP::Atom *Md2QmEngine::AddAtomType(CTP::Molecule *owner, int residue_number,
                                    string md_atom_name,  int md_atom_id,
                                    bool hasQMPart,       int qm_atom_id,
                                    double weight) {
    CTP::Atom* atom = new CTP::Atom(owner,                residue_number,
                                    md_atom_name,         md_atom_id,
                                    hasQMPart,            qm_atom_id,
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

CTP::Segment *Md2QmEngine::AddSegmentType(int segment_id, Property *property) {
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


CTP::Molecule *Md2QmEngine::getMoleculeType(const string &name) {
    try {
        return _map_MoleculeName_MoleculeType.at(name);
    }
    catch ( out_of_range ) {
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
    catch ( out_of_range ) {
        return mdname;
    }
}

CTP::Atom *Md2QmEngine::getAtomType(const string &molMdName,
                                    int resNr, const string &mdAtomName) {
    return this->_map_mol_resNr_atm_atmType.at(molMdName)
                                           .at(resNr)
                                           .at(mdAtomName);
}


void Md2QmEngine::PrintInfo() {

    vector<CTP::Molecule*>::iterator mit;
    vector<CTP::Segment*>::iterator sit;
    vector<CTP::Fragment*>::iterator fit;
    vector<CTP::Atom*>::iterator ait;

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

    for (mit = _molecule_types.begin();
         mit != _molecule_types.end();
         mit++) {
         CTP::Molecule *mol = *mit;
         cout << "[ " << mol->getName() << " ]" << endl;

         for (sit = mol->Segments().begin();
              sit != mol->Segments().end();
              sit++) {
              CTP::Segment *seg = *sit;
              cout << " - Segment [ " << seg->getName() << " ]" << endl;

              for (fit = seg->Fragments().begin();
                   fit != seg->Fragments().end();
                   fit++ ) {
                   CTP::Fragment *frag = *fit;
                   cout << "   - Fragment [ " << frag->getName() << " ]: "
                      << frag->Atoms().size() << " atoms " << endl;
              }
         }
    }

    map < string, string > ::iterator mssit;
    for (mssit = this->_map_MoleculeMDName_MoleculeName.begin();
         mssit != this->_map_MoleculeMDName_MoleculeName.end();
         mssit ++) {
         cout << "MD [ " << mssit->first << " ] mapped to "
              << "QM [ " << mssit->second << " ] " << endl;
    }

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
                  cout << "MD mol.name " << it0->first
                       << " | MD res.number " << it1->first
                       << " | MD atm.name " << it2->first
                       << " => QM atm.id " << it2->second->getQMId()
                       << endl;
             }
         }
    }


}



