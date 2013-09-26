#ifndef VOTCA_CTP_PEWALD3D_H
#define VOTCA_CTP_PEWALD3D_H

#include <votca/csg/boundarycondition.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/ewaldnd.h>

namespace CSG = votca::csg;

namespace votca { namespace ctp {
    
    // NOTE: This is not a conventional 3D Ewald summation, so use carefully
    //       (tuned for the purpose of site-energy calculations)
    // NOTE: PolarTop should be set-up with three containers: FGC, FGN, BGN
    //       MGN is set-up in constructor using the real-space c/o (from input)
    //       The topology is used to retrieve information on the sim. box (PB)
    //       All polar segments should be positioned as nearest images of the
    //       foreground charge density (FGC, FGN).
    //       All polar segments should be appropriately charged (Q00).
    // NOTE: The k-shell grouping algorithm can fail for strongly skewed boxes.
    //        
    
    class PEwald3D3D : public Ewald3DnD
    {
        
    public:
        
        PEwald3D3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
       ~PEwald3D3D();
        
        string IdentifyMethod() { return "Polar 3D x 3D"; }
        void GenerateKVectors();
        
        EWD::triple<> ConvergeRealSpaceSum();
        EWD::triple<> ConvergeReciprocalSpaceSum();
        EWD::triple<> CalculateForegroundCorrection();
        EWD::triple<> CalculateShapeCorrection();
        EWD::triple<> CalculateHigherRankCorrection() { return EWD::triple<>(0,0,0); }        
        EWD::triple<> CalculateK0Correction() { return EWD::triple<>(0,0,0); }
        
        void Field_ConvergeRealSpaceSum();
        void Field_ConvergeReciprocalSpaceSum();
        void Field_CalculateForegroundCorrection();
        void Field_CalculateShapeCorrection();
    
    private:
        
        string _shape; // Summation shape (for 3D corr. term)
        
    };

}}


#endif