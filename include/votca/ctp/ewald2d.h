#ifndef VOTCA_CTP_EWALD2D_H
#define VOTCA_CTP_EWALD2D_H

#include <votca/csg/boundarycondition.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/xinteractor.h>

namespace CSG = votca::csg;

namespace votca { namespace ctp {
    
    // NOTE: This is not a conventional 2D Ewald summation, so use carefully
    //       (tuned for the purpose of site-energy calculations)
    
    class Ewald2D
    {
        
    public:
        
        Ewald2D(Topology *top, PolarTop *ptop, XJob *xjob, Logger *log);
       ~Ewald2D() { ; }
       
        void Evaluate();
        
        
        
    private:
        
        XInteractor _actor;
        Logger *_log;
        
        // PERIODIC BOUNDARY
        Topology *_top;
        CSG::BoundaryCondition *_bc;    // Periodicity reduced to xy sub-space
        
        // POLAR SEGMENTS
        PolarTop *_ptop;
        vector< PolarSeg* > _bg_P;      // Period. density = _bg_N v _fg_N
        vector< PolarSeg* > _bg_N;      // Neutral background
        vector< PolarSeg* > _mg_N;      // Neutral midground
        vector< PolarSeg* > _fg_N;      // Neutral foreground
        vector< PolarSeg* > _fg_C;      // Charged foreground
        
        // CONVERGENCE
        double _a;                      // _a = 1/(sqrt(2)*sigma)
        double _K_co;                   // k-space c/o
        double _R_co;                   // r-space c/o
        
    };    
    
}}


#endif