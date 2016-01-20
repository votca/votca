#ifndef VOTCA_XTP_POISSONGRID_H
#define VOTCA_XTP_POISSONGRID_H

#include <votca/xtp/polartop.h>
#include <votca/xtp/logger.h>

namespace votca { 
namespace xtp {
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