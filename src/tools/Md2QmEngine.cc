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

    // Load typology from xml-file to fill
    // mol./segm./fragm./atom type containers

    Property typology;
    load_property_from_xml(typology, xmlfile.c_str());

    string key = "topology.molecules.molecule";
    list<Property *> molecules = typology.Select(key);
    list<Property *>::iterator it_molecule;

    // cout << " Found " << molecules.size() << " molecule type(s)" << endl;

    // +++++++++++++++++++++++++++++++++++++++++++ //
    // Molecules => Segments => Fragments => Atoms //
    // +++++++++++++++++++++++++++++++++++++++++++ //

    // Load molecules, create fragments and segments
    int molecule_id = 1;
    for ( it_molecule = molecules.begin();
           it_molecule != molecules.end();
           ++it_molecule ) {

       CTP::Molecule *molecule = AddMoleculeType(molecule_id++, *it_molecule );

       // Load segments
       key = "segments.segment";
       list<Property *> segments = (*it_molecule)->Select(key);
       list<Property *>::iterator it_segment;
       //cout << "  - Found " << segments.size()
       //     << " segment(s) in molecule type " << molecule->getName() << endl;

       int segment_id = 1;

       for ( it_segment = segments.begin();
              it_segment != segments.end();
              ++it_segment ) {

         // Create a new segment
         CTP::Segment *segment = AddSegmentType(segment_id++, *it_segment );
         molecule->AddSegment( segment );

         // Load fragments
         key = "fragments.fragment";
         list<Property *> fragments = (*it_segment)->Select(key);
         list<Property *>::iterator it_fragment;
         //cout << "    - Found " << fragments.size()
         //         << " fragments in segment " << segment->getName() << endl;

          int fragment_id = 1;
          for ( it_fragment = fragments.begin();
                 it_fragment != fragments.end();
                 ++it_fragment ) {

             CTP::Fragment* fragment = AddFragmentType(fragment_id,
                                                      *it_fragment);
             segment->AddFragment( fragment );

             string mdatoms = (*it_fragment)->get("mdatoms").as<string>();
             string qmatoms = (*it_fragment)->get("qmatoms").as<string>();
             string weights = (*it_fragment)->get("weights").as<string>();

             Tokenizer tok_md_atoms(mdatoms, " ");
             Tokenizer tok_qm_atoms(qmatoms, " ");
             Tokenizer tok_weights(weights, " ");

             vector <string> md_atom_names;
             vector <string> qm_atom_names;
             vector <string> atom_weights;;

             tok_md_atoms.ToVector(md_atom_names);
             tok_weights.ToVector(atom_weights);

             vector<string>::iterator it_md_atom_name;
             vector<string>::iterator it_atom_weight;
             for ( it_md_atom_name = md_atom_names.begin(),
                   it_atom_weight  = atom_weights.begin();
                    it_md_atom_name != md_atom_names.end();
                    ++it_md_atom_name, ++ it_atom_weight) {

               // Process atom info
               Tokenizer tok_md((*it_md_atom_name), ":");
               vector<string> md_atom_info;
               tok_md.ToVector( md_atom_info );
               int residue_number = boost::lexical_cast<int>(md_atom_info[0]);
               string residue_name = md_atom_info[1];
               string md_atom_name =  md_atom_info[2];
               double weight = boost::lexical_cast<double>(*it_atom_weight);

               // TODO: QM atoms + atom_ids
               int atom_id = 0;

               CTP::Atom *atom = AddAtomType(molecule, atom_id, md_atom_name,
                                             residue_number, weight);
               fragment->AddAtom(atom);
               segment->AddAtom(atom);
               molecule->AddAtom(atom);

             }
             //cout << "      - Found " << fragment->NumberOfAtoms()
             //     << " atoms in fragment " <<  fragment->getName() << endl;
          }
          //cout <<  "    - Total of " << segment->NumberOfAtoms()
          //     << " atoms in segment " << segment->getName() << endl;

       }
       //cout <<  "  - Total of " << molecule->NumberOfAtoms()
       //     << " atoms in the molecule " << molecule->getName() << endl;
    }

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

        CSG::Molecule *molMD = *mit;
        ++globalMolID;
        // Use molecule name to create associated QM|MD molecule, segments, ...
        // with yet empty atom containers

        for (int i = 0; i < molMD->BeadCount(); i++) {
            CSG::Bead* atomMD = molMD->getBead(i);








            cout << "ATOM getId  " << atomMD->getId();
            cout << " | getName  " << atomMD->getName();
            cout << " | getMol   " << atomMD->getMolecule()->getId();
            cout << " | MolName  " << atomMD->getMolecule()->getName();
            cout << " | getResnr " << atomMD->getResnr();
            cout << " | getType  " << atomMD->getType()->getName();
            cout << " | getPos   " << atomMD->getPos();
            cout << endl;
        }
    }
}















CTP::Atom *Md2QmEngine::AddAtomType(CTP::Molecule *owner, int atom_id,
                                    string atom_name, int residue_number,
                                    double weight) {
    CTP::Atom* atom = new CTP::Atom(owner, atom_id, atom_name,
                                    residue_number, weight);
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
    _map_MoleculeMDName_MoleculeName[molecule_name] = molecule_mdname;
    _map_MoleculeName_MoleculeType[molecule_name] = molecule;
    return molecule;
}


CTP::Molecule *Md2QmEngine::getMoleculeType(const string &name) {
    return _map_MoleculeName_MoleculeType[name];
}

const string &Md2QmEngine::getMoleculeName(const string &mdname) {
    return _map_MoleculeMDName_MoleculeName[mdname];
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


}



