#ifndef __POLARSEG__H
#define	__POLARSEG__H

#include <votca/tools/vec.h>
#include <votca/ctp/apolarsite.h>

namespace votca { namespace ctp {
    
class PolarSeg : public vector<APolarSite*>
{

public:

    PolarSeg(int id, vector<APolarSite*> &psites);
   ~PolarSeg();

    const int getId() { return _id; }
    const vec &getPos() { return _pos; }
    void CalcPos();
    double CalcTotQ();
    void Translate(const vec &shift);

private:

    int _id;
    vec _pos;    


};    
    
    
}}

#endif