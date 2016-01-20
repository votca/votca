#ifndef VOTCA_XTP_MOLPOLENGINE_H
#define VOTCA_XTP_MOLPOLENGINE_H

#include <votca/xtp/polartop.h>

namespace votca {
namespace xtp {
    
    
class MolPolEngine
{
public:
    
    MolPolEngine() 
        : _aDamp(0.390), _wSOR(0.30), _maxIter(1024), _epsTol(0.0001)
        { _actor.SetADamp(_aDamp); }
    MolPolEngine(double aDamp, double wSOR, int maxIter, double epsTol)
        : _aDamp(aDamp), _wSOR(wSOR), _maxIter(maxIter), _epsTol(epsTol)
        { _actor.SetADamp(_aDamp); }
   ~MolPolEngine() {}
    
    matrix CalculateMolPol(vector<APolarSite*> &poles, bool verbose = true);
    int    SCF_Induce(vector<APolarSite*> &poles);
    
private:
    
    BasicInteractor _actor;
    double _aDamp;
    double _wSOR;
    int _maxIter;
    double _epsTol;
};

}}

#endif