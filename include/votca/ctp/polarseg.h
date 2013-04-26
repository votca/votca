#ifndef __POLARSEG__H
#define	__POLARSEG__H

#include <votca/tools/vec.h>
#include <votca/ctp/apolarsite.h>

namespace votca { namespace ctp {
    
class PolarSeg : public vector<APolarSite*>
{

public:

    PolarSeg(vector<APolarSite*> &psites);
   ~PolarSeg();

    const vec &getPos() { return _pos; }
    void CalcPos();
    void Translate(const vec &shift);

private:

    vec _pos;    


};    
    
    
}}

#endif