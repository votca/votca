#ifndef VOTCA_CTP_POLARFRAG_H
#define VOTCA_CTP_POLARFRAG_H

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
    const vec &CalcPos() {    
        _pos = vec(0,0,0);    
        for (int i = 0; i < this->size(); ++i) _pos += (*this)[i]->getPos();
        if (this->size() > 0) _pos /= double(this->size());
        return _pos;
    }
    
    
private:
    int _id;
    PolarSeg *_pseg;
    vec _pos;
};

}}

#endif