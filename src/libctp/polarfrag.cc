#include <votca/ctp/polarfrag.h>

namespace votca {
namespace ctp {


const votca::tools::vec &PolarFrag::CalcPosCenterOfGeom() {    
    _pos = votca::tools::vec(0,0,0);    
    for (int i = 0; i < this->size(); ++i) {
        _pos += (*this)[i]->getPos();
    }
    if (this->size() > 0) _pos /= double(this->size());
    return _pos;
}
    
    
const votca::tools::vec &PolarFrag::CalcPosPolarWeights() {
    // Establish total weights
    double sum_iso_p = 0.0;
    double sum_abs_q = 0.0;
    for (int i = 0; i < this->size(); ++i) {
        sum_iso_p += (*this)[i]->getIsoP();
        sum_abs_q += std::abs((*this)[i]->getQ00());
    }

    _pos = votca::tools::vec(0,0,0);
    // Any polar sites in here? If not, return (0,0,0)
    if (this->size() < 1) ;
    // Noteworthy charge and polarizability? If not, return CoG
    else if (sum_abs_q <= 1e-2 && sum_iso_p <= 1e-2) {
        _pos = this->CalcPosCenterOfGeom();
    }
    // Else: return weighted 0.5*(center-of-polarity + center-of-charge)
    else {
        for (int i = 0; i < this->size(); ++i) {
            double weight = 0.0;
            if (sum_abs_q < 1e-2) {
                assert(sum_iso_p >= 1e-2 && "<CalcPosPolarWeights> P-ERROR");
                weight = (*this)[i]->getIsoP()/sum_iso_p;
            }
            else if (sum_iso_p < 1e-2) {
                assert(sum_abs_q >= 1e-2 && "<CalcPosPolarWeights> Q-ERROR");
                weight = std::abs((*this)[i]->getQ00())/sum_abs_q;
            }
            else {
                weight = 0.5 * (
                    (*this)[i]->getIsoP()/sum_iso_p
                  + std::abs((*this)[i]->getQ00())/sum_abs_q);
            }
            _pos += (*this)[i]->getPos() * weight;
        }
    }
    return _pos;
}
    
    
    
    
    
    
    
    
    
    
    
}}