#ifndef __POLARSEG__H
#define	__POLARSEG__H

#include <votca/tools/vec.h>
#include <votca/ctp/apolarsite.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace votca { namespace ctp {
    
class PolarSeg : public vector<APolarSite*>
{

public:

    PolarSeg() {};
    PolarSeg(int id, vector<APolarSite*> &psites);
    PolarSeg(PolarSeg *templ);
   ~PolarSeg();

    const int &getId() { return _id; }
    const vec &getPos() { return _pos; }
    void       CalcPos();
    double     CalcTotQ();
    void       Translate(const vec &shift);
    
    template<class Archive>
    void serialize(Archive &arch, const unsigned int version) {
        arch & boost::serialization::base_object< vector<APolarSite*> >(*this);
        arch & _id;
        arch & _pos;
        return;
    }
    
private:

    int _id;
    vec _pos;    


};


}}

#endif