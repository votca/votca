#include <votca/ctp/polarseg.h>


namespace votca { namespace ctp {

    
    
PolarSeg::PolarSeg(int id, vector<APolarSite*> &psites) : _id(id) {    
    for (int i = 0; i < psites.size(); ++i) {
        push_back(psites[i]);
    }    
    this->CalcPos();
}
    
    
PolarSeg::~PolarSeg() {
   vector<APolarSite*> ::iterator pit;
   for (pit = begin(); pit < end(); ++pit) {         
       delete *pit;
   }
   clear();
}


void PolarSeg::CalcPos() {    
    _pos = vec(0,0,0);    
    for (int i = 0; i < this->size(); ++i) {        
        _pos += (*this)[i]->getPos();        
    }
    _pos /= double(this->size());      
}


double PolarSeg::CalcTotQ() {
    double Q = 0.0;
    for (int i = 0; i < this->size(); ++i) {
        Q += (*this)[i]->getQ00();
    }
    return Q;
}


void PolarSeg::Translate(const vec &shift) {    
    for (int i = 0; i < size(); ++i) {
        (*this)[i]->Translate(shift);
    }
    _pos += shift;
}



}}