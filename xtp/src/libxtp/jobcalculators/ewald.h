#ifndef VOTCA_CTP_EWALD_H
#define VOTCA_CTP_EWALD_H


#include <votca/xtp/parallelxjobcalc.h>
//#include <votca/ctp/xmapper.h>
#include <votca/xtp/ewald/xjob.h>
#include <votca/xtp/ewald/xinductor.h>
#include <votca/xtp/ewald/xinteractor.h>
//#include <votca/ctp/ewald2d.h>
//#include <votca/ctp/ewald3d.h>
#include <votca/xtp/ewald/pewald3d.h>
#include <votca/xtp/logger.h>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>


using boost::format;


namespace votca { namespace xtp {


template<class EwaldMethod>
class Ewald : public ParallelXJobCalc< std::vector<Job>>
{
public:
    
    Ewald() {};
   ~Ewald() {};
   
    std::string          Identify() const { return "ewald"; }
    void            Initialize(tools::Property *);
    void            WriteJobFile(const Topology& top);
    void            ReadJobFile(Topology& top);
    
    void            PreProcess(Topology *top);
    Job::JobResult  EvalJob(const Topology& top, Job& job, QMThread& Thread);
    void            PostProcess(Topology *top) { ; }
    
    XJob            ProcessInputString(Job &, Topology &, QMThread &);   

 protected:
  void ParseSpecificOptions(const tools::Property& user_options); 

private:
    
    tools::Property                      *_options;
    // MULTIPOLES DEFINITION & MAPPING
    std::string                         _xml_file;
    std::string                         _mps_table;    
    std::string                         _polar_bg_arch;
    //XMpsMap                        _mps_mapper;
    bool                           _pdb_check;
    bool                           _ptop_check;
};




}}

#endif

