#ifndef VOTCA_CTP_EWALDND_H
#define VOTCA_CTP_EWALDND_H

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
    
    class Ewald3DnD
    {
        
    public:
        
        Ewald3DnD() { ; }
        Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
        virtual ~Ewald3DnD();
       
        void SetupMidground(double R_co);
        void WriteDensitiesPDB(string pdbfile);
        void Evaluate();
        bool Converged() { return _converged_R && _converged_K; }
        string GenerateOutputString();
        string GenerateErrorString();
        
        virtual string IdentifyMethod() = 0;
        virtual double ConvergeRealSpaceSum();
        virtual double ConvergeReciprocalSpaceSum() = 0;
        virtual double CalculateForegroundCorrection();
        virtual double CalculateHigherRankCorrection();
        virtual double CalculateShapeCorrection() { return 0.0; }
        virtual double CalculateK0Correction() { return 0.0; }
        
        // To sort K-vectors via std::sort using a norm functor
        template<class Norm>
        struct VectorSort
        {
            VectorSort() : _p(1e-40) { ; }
            VectorSort(double precision) : _p(precision) { ; }
            inline bool operator() (const vec &v1, const vec &v2);
            inline bool MatchDouble(double a, double b) 
                { return ((a-b)*(a-b) < _p) ? true : false; }
            double _p;
            Norm _norm;
        };
        
        // Tschebyschow norm functor
        struct MaxNorm { inline double operator() (const vec &v) 
            { return votca::tools::maxnorm(v); } };
        // Euclidean norm functor
        struct EucNorm { inline double operator() (const vec &v) 
            { return votca::tools::abs(v); } };
        
    protected:
        
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
        
        VectorSort<MaxNorm> _maxsort;
        VectorSort<EucNorm> _eucsort;
        
        // ENERGIES
        double _ER;                     // R-space sum
        double _EC;                     // R-space correction
        double _EK;                     // K-space sum
        double _E0;                     // K-space K=0 contribution
        double _ET;                     // ER - EC + EK + E0
        double _EDQ;                    // Higher-Rank FGC->MGN correction
        double _EJ;                     // Geometry-dependent correction
        
        
        
    };


template<class Norm>
inline bool Ewald3DnD::VectorSort<Norm>::operator() (const vec &v1,
    const vec &v2) {
    bool smaller = false;
    // LEVEL 1: MAGNITUDE
    double V1 = _norm(v1);
    double V2 = _norm(v2);
    if (MatchDouble(V1,V2)) {
        // LEVEL 2: X
        double X1 = v1.getX();
        double X2 = v2.getX();
        if (MatchDouble(X1,X2)) {
            // LEVEL 3: Y
            double Y1 = v1.getY();
            double Y2 = v2.getY();
            if (MatchDouble(Y1,Y2)) {
                // LEVEL 4: Z
                double Z1 = v1.getZ();
                double Z2 = v2.getZ();
                if (MatchDouble(Z1,Z2)) smaller = true;
                else smaller = (Z1 < Z2) ? true : false;
            }
            else smaller = (Y1 < Y2) ? true : false;
        }
        else smaller = (X1 < X2) ? true : false;
    }
    else smaller = (V1 < V2) ? true : false;          
    return smaller;
}


}}


#endif