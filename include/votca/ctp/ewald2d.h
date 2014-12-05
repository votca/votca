#ifndef VOTCA_CTP_EWALD2D_H
#define VOTCA_CTP_EWALD2D_H

#include <votca/csg/boundarycondition.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/ewaldnd.h>

namespace CSG = votca::csg;

namespace votca { namespace ctp {
    
    // NOTE: This is not a conventional 2D Ewald summation, so use carefully
    //       (tuned for the purpose of site-energy calculations)
    // NOTE: PolarTop should be set-up with three containers: FGC, FGN, BGN
    //       MGN is set-up in constructor using the real-space c/o (from input)
    //       The topology is used to retrieve information on the sim. box (PB)
    //       All polar segments should be positioned as nearest images of the
    //       foreground charge density (FGC, FGN).
    
    class Ewald3D2D : public Ewald3DnD
    {
        
    public:
        
        Ewald3D2D(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
       ~Ewald3D2D();
       
        string IdentifyMethod() { return "3D x 2D"; }
        EWD::triple<> ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target);
        EWD::triple<> CalculateK0Correction(vector<PolarSeg*> &target);
                
    };
    
}}


#endif