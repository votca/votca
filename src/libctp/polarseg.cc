#include <votca/ctp/polarseg.h>


namespace votca { namespace ctp {

    
    
PolarSeg::PolarSeg(vector<APolarSite*> &psites) : _psites(psites) {    
    this->CalcPos();
}
    
    
PolarSeg::~PolarSeg() {
   vector<APolarSite*> ::iterator pit;
   for (pit = _psites.begin(); pit < _psites.end(); ++pit) {               
       delete *pit;               
   }
   _psites.clear();
}


void PolarSeg::CalcPos() {    
    _pos = vec(0,0,0);    
    for (int i = 0; i < _psites.size(); ++i) {        
        _pos += _psites[i]->getPos();        
    }    
    _pos /= double(_psites.size());      
}


void PolarSeg::Translate(const vec &shift) {    
    for (int i = 0; i < _psites.size(); ++i) {
        _psites[i]->Translate(shift);
    }
    _pos += shift;
}



}}