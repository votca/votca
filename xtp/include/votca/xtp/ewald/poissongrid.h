#ifndef VOTCA_CTP_POISSONGRID_H
#define VOTCA_CTP_POISSONGRID_H

#include <votca/xtp/ewald/polartop.h>
#include <votca/xtp/logger.h>

namespace votca { 
namespace xtp {
namespace POI {
    
   
class PoissonGrid
{
public:
    
    // Also take into account shape
    PoissonGrid(const Topology *top, std::vector<PolarSeg*> &fg, std::vector<PolarSeg*> &bg, Logger* log);
    
    
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