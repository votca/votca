#ifndef VOTCA_CTP_EWALD3D_H
#define VOTCA_CTP_EWALD3D_H

#include <votca/csg/boundarycondition.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/xinteractor.h>

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
    
    class Ewald3D
    {
        
    public:
        
        Ewald3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
       ~Ewald3D();
       
        void SetupMidground(double R_co);
        void WriteDensitiesPDB(string pdbfile);
        
        void Evaluate();
        double ConvergeRealSpaceSum();
        double ConvergeReciprocalSpaceSum();
        double CalculateShapeCorrection(string shape);
        double CalculateSq2(vec &k);
        
        bool Converged() { return _converged_R && _converged_K; }
        string GenerateOutputString();
        string GenerateErrorString();
        
        // To sort K-vectors via std::sort - Euclidean norm
        struct EuclideanSort
        {
            inline bool operator() (const vec &v1, const vec &v2, double t=1e-40);            
            inline bool MatchDouble(double a, double b, double t=1e-40);
        };
        
        // To sort K-vectors via std::sort - Tschebyschow norm
        struct TschebyschowSort
        {
            inline bool operator() (const vec &v1, const vec &v2, double t=1e-40);
            inline bool MatchDouble(double a, double b, double t=1e-40);
        };
        
        
    private:
        
        XInteractor _actor;
        Logger *_log;
        
        // PERIODIC BOUNDARY
        Topology *_top;
        CSG::BoundaryCondition *_bc;    // Periodicity reduced to xy sub-space
        vec _center;
        
        // POLAR SEGMENTS
        PolarTop *_ptop;
        vector< PolarSeg* > _bg_P;      // Period. density = _bg_N v _fg_N
        vector< PolarSeg* > _bg_N;      // Neutral background
        vector< PolarSeg* > _mg_N;      // Neutral midground
        vector< PolarSeg* > _fg_N;      // Neutral foreground
        vector< PolarSeg* > _fg_C;      // Charged foreground
        vector< bool > _inForeground;
        
        // CONVERGENCE
        double _alpha;                  // _a = 1/(sqrt(2)*sigma)
        double _K_co;                   // k-space c/o
        double _R_co;                   // r-space c/o
        double _crit_dE;                // Energy convergence criterion [eV]
        bool   _converged_R;            // Did R-space sum converge?
        bool   _converged_K;            // Did K-space sum converge?
        
        // LATTICE (REAL, RECIPROCAL)
        vec _a; vec _b; vec _c;         // Real-space lattice vectors
        int _na_max, _nb_max, _nc_max;  // Max. cell indices to sum over (R)
        vec _A; vec _B; vec _C;         // Reciprocal-space lattice vectors
        int _NA_max, _NB_max, _NC_max;  // Max. cell indices to sum over (K)
        double _LxLy;                   // |a^b|
        double _LxLyLz;                 // a*|b^c|
        string _shape;                  // Summation shape (for 3D corr. term)
        
        TschebyschowSort _maxsort;
        EuclideanSort _eucsort;
        
        // ENERGIES
        double _ER;                     // R-space sum
        double _EC;                     // R-space correction
        double _EK;                     // K-space sum
        double _E0;                     // K-space K=0 contribution
        double _ET;                     // ER - EC + EK + E0
        double _EDQ;                    // Higher-Rank FGC->MGN correction
        double _EJ;                     // Geometry-dependent correction
        
        
        
    };


inline bool Ewald3D::EuclideanSort::operator() (const vec& v1, 
    const vec& v2, double t) {
    bool smaller = false;
    // LEVEL 1: MAGNITUDE
    double V1 = abs(v1);
    double V2 = abs(v2);                
    if (MatchDouble(V1,V2,t)) {
        // LEVEL 2: X
        double X1 = v1.getX();
        double X2 = v2.getX();
        if (MatchDouble(X1,X2,t)) {
            // LEVEL 3: Y
            double Y1 = v1.getY();
            double Y2 = v2.getY();
            if (MatchDouble(Y1,Y2,t)) {
                // LEVEL 4: Z
                double Z1 = v1.getZ();
                double Z2 = v2.getZ();
                if (MatchDouble(Z1,Z2,t)) smaller = true;
                else smaller = (Z1 < Z2) ? true : false;
            }
            else smaller = (Y1 < Y2) ? true : false;
        }
        else smaller = (X1 < X2) ? true : false;
    }
    else smaller = (V1 < V2) ? true : false;          
    return smaller;
}


inline bool Ewald3D::EuclideanSort::MatchDouble(double a, double b, 
    double t) {
    return ((a-b)*(a-b) < t) ? true : false;
}


inline bool Ewald3D::TschebyschowSort::operator() (const vec& v1, 
    const vec& v2, double t) {
    bool smaller = false;
    // LEVEL 1: MAGNITUDE
    double V1 = maxnorm(v1);
    double V2 = maxnorm(v2);
    if (MatchDouble(V1,V2,t)) {
        // LEVEL 2: X
        double X1 = v1.getX();
        double X2 = v2.getX();
        if (MatchDouble(X1,X2,t)) {
            // LEVEL 3: Y
            double Y1 = v1.getY();
            double Y2 = v2.getY();
            if (MatchDouble(Y1,Y2,t)) {
                // LEVEL 4: Z
                double Z1 = v1.getZ();
                double Z2 = v2.getZ();
                if (MatchDouble(Z1,Z2,t)) smaller = true;
                else smaller = (Z1 < Z2) ? true : false;
            }
            else smaller = (Y1 < Y2) ? true : false;
        }
        else smaller = (X1 < X2) ? true : false;
    }
    else smaller = (V1 < V2) ? true : false;          
    return smaller;
}

inline bool Ewald3D::TschebyschowSort::MatchDouble(double a, double b, 
    double t) {
    return ((a-b)*(a-b) < t) ? true : false;
}

}}


#endif