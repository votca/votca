#include <votca/ctp/fragment.h>

using namespace std;
namespace votca { namespace ctp {

Fragment::~Fragment() {
    vector < Atom* > ::iterator atmit;
    for (atmit = this->Atoms().begin();
          atmit < this->Atoms().end();
          atmit++) {
          delete *atmit;
    }
    _weights.clear();
    _atoms.clear();
}


void Fragment::AddAtom(Atom* atom) {
    _atoms.push_back( atom );
    atom->setFragment(this);
    _weights.push_back( atom->getWeight() );
}

void Fragment::Rotate() {
    throw runtime_error("Not implemented.");
}

void Fragment::Translate() {
    throw runtime_error("Not implemented.");
}

void Fragment::calcPos() {
    vec pos = vec(0,0,0);
    double totWeight = 0.0;

    for (int i = 0; i< _atoms.size(); i++) {
        pos += _atoms[i]->getPos() * _atoms[i]->getWeight();
        totWeight += _atoms[i]->getWeight();
    }

    pos = pos / totWeight;

}


}}