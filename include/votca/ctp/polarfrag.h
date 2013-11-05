#ifndef VOTCA_CTP_POLARFRAG_H
#define VOTCA_CTP_POLARFRAG_H

#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace votca { namespace ctp {
    
    
// Forward declarations
class PolarSeg;
class APolarSite;

class PolarFrag : public std::vector<APolarSite*>
{    
public:    
    PolarFrag(PolarSeg *pseg, int id) : _id(id), _pseg(pseg) {}
   ~PolarFrag() { /* Polar sites cleaned by PolarSeg */; }
    int getId() { return _id; }
    const vec &CalcPosCenterOfGeom() {    
        _pos = vec(0,0,0);    
        for (int i = 0; i < this->size(); ++i) {
            _pos += (*this)[i]->getPos();
        }
        if (this->size() > 0) _pos /= double(this->size());
        return _pos;
    }
    const vec &CalcPosPolarWeights() {    
        
        double sum_iso_p = 0.0;
        double sum_abs_q = 0.0;
        for (int i = 0; i < this->size(); ++i) {
            sum_iso_p += (*this)[i]->getIsoP();
            sum_abs_q += std::abs((*this)[i]->getQ00());
        }
        
        _pos = vec(0,0,0);
        // Any polar sites in here? If not, return (0,0,0)
        if (this->size() < 1)
            ;
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
    
    // Serialization interface
    template<class Archive>
    void serialize(Archive &arch, const unsigned int version) {
        arch & boost::serialization::base_object< vector<APolarSite*> >(*this);
        arch & _id;
        arch & _pseg;
        arch & _pos;
        return;
    }
    
    
private:
    int _id;
    PolarSeg *_pseg;
    vec _pos;
};

}}

#endif