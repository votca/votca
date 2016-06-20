
#include <votca/xtp/nbo.h>
#include <votca/xtp/aomatrix.h>
#include <votca/tools/linalg.h>
#include <boost/progress.hpp>
#include <votca/xtp/numerical_integrations.h>
#include <math.h> 


using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    

void Nbo::Evaluate(std::vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs) { 
   
    //TODO: Jens, fill this in
    
    return;
    }


void Nbo::LoadMatrices(std::string fn_projectionMatrix, std::string fn_overlapMatrix){
    
    //TODO: Yuriy, fill this in
    
    return;
    }

}}