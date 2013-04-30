#include <votca/ctp/qmmachine.h>


namespace votca { namespace ctp {


QMMachine::QMMachine(XJob *job, XInductor *xind, Gaussian *qmpack,
                     Property *opt, string sfx, bool mav) 
                   : _job(job), _xind(xind), _qmpack(qmpack) {
    
    string key = sfx + ".qmmm";    
    
}
   
    
    
    
    
    
}}
