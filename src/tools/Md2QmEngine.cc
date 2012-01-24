#include "Md2QmEngine.h"

using namespace votca::ctp;
using namespace votca::tools;
/*
Md2QmEngine::~Md2QmEngine() {

    vector < Fragment* > :: iterator fragment;
    vector < Segment* > :: iterator segment;
    vector < Molecule* > :: iterator molecule; 

    // Clean up list of fragment types
    cout << " * Typology: deleting fragment types: ";
    for (fragment = _fragment_types.begin();
         fragment < _fragment_types.end();
         ++fragment) {
         cout << (*fragment)->getId() << ":"
              << (*fragment)->getName() << " ";
         delete *fragment;
    }
    cout << endl;
    _fragment_types.clear();

    // Clean up list of segment types
    cout << " * Typology: deleting segment types: ";
    for (segment = _segment_types.begin();
         segment < _segment_types.end();
         ++segment) {
         cout << (*segment)->getId() << ":"
              << (*segment)->getName() << " ";
         delete *segment;
    }
    cout << endl;
    _segment_types.clear();

    // Clean up list of molecule types; this also deletes atoms
    cout << " * Typology: deleting molecule types: ";
    for (molecule = _molecule_types.begin();
         molecule < _molecule_types.end(); ++molecule){
         cout << (*molecule)->getId() << ":"
              << (*molecule)->getName() << " ";
         delete *molecule;
    }
    cout << endl;
    _molecule_types.clear();

    // clean up the map of the molecule md name : pointer to molecule type
    _map_MoleculeName_MoleculeType.clear();
    _map_MoleculeMDName_MoleculeName.clear();

}
*/

void Md2QmEngine::Initialize(const string &xmlfile) {

    cout << "Md2QmEngine::Initialize ..." << endl;
    cout << "Mapping taken from " << xmlfile << endl;

}
