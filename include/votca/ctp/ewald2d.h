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
        
        Ewald2D(Topology *top, PolarTop *ptop, XJob *xjob, double R_co, Logger *log);
       ~Ewald2D() { ; }
       
        void CheckParameters();
        void Evaluate();
        //bool SortByMagnitude(vec k1, vec k2) { return abs(k1) < abs(k2); }
        
        
        struct vec_smaller_than
        {
            inline bool operator() (const vec& vec1, const vec& vec2, double t = 1e-40)
            {
                bool smaller = false;                
                // LEVEL 1: Magnitude
                double V1 = abs(vec1);
                double V2 = abs(vec2);                
                if (match_double(V1,V2,t)) {
                    // LEVEL 2: X
                    double X1 = vec1.getX();
                    double X2 = vec2.getX();
                    if (match_double(X1,X2,t)) {
                        // LEVEL 3: Y
                        double Y1 = vec1.getY();
                        double Y2 = vec2.getY();
                        if (match_double(Y1,Y2,t)) {
                            // LEVEL 4: Z
                            double Z1 = vec1.getZ();
                            double Z2 = vec2.getZ();
                            if (match_double(Z1,Z2,t)) smaller = true;
                            else if (Z1 < Z2) smaller = true;
                            else smaller = false;
                        }
                        else if (Y1 < Y2) smaller = true;
                        else smaller = false;
                    }
                    else if (X1 < X2) smaller = true;
                    else smaller = false;
                }
                else if (V1 < V2) smaller = true;
                else smaller = false;                
                return smaller;
            }
            
            inline bool match_double(double a, double b, double t = 1e-40) {                
                if ((a-b)*(a-b) < t) {
                    return true;
                }
                else {
                    return false;
                }
            }
            
        };
        
        
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