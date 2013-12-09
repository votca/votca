#ifndef VOTCA_CTP_POISSONGRID_H
#define VOTCA_CTP_POISSONGRID_H

#include <votca/ctp/polartop.h>
#include <votca/ctp/logger.h>

namespace votca { 
namespace ctp {
namespace POI {
    
   
class PoissonGrid
{
public:
    
    // Also take into account shape
    PoissonGrid(Topology *top, vector<PolarSeg*> &fg, vector<PolarSeg*> &bg, Logger* log);
    
    
private:
    
    Logger *_log;
    vec _center;
    
    
};



class PoissonCell
{
    
    PoissonCell();
    
    
    
};
    
    
    
    
    
    
    
    
    
}}}

#endif