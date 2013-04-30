#ifndef __QMMACHINE__H
#define	__QMMACHINE__H


#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/gaussian.h>


namespace votca { namespace ctp {
    
class QMMachine 
{

public:

    QMMachine(XJob *job, XInductor *xind, Gaussian *qmpack,
              Property *opt, string sfx, bool mav);
   ~QMMachine() {};   
   
    void Evaluate(XJob *job);

private:

    XJob *_job;
    XInductor *_xind;
    Gaussian *_qmpack;

};

    
}}

#endif