#ifndef VOTCA_CTP_EWALDND_H
#define VOTCA_CTP_EWALDND_H

#include <votca/csg/boundarycondition.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/xinteractor.h>
#include <votca/ctp/ewaldactor.h>
#include <votca/ctp/xinductor.h>

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
       
        void ExpandForegroundReduceBackground(double polar_R_co);
        void SetupMidground(double R_co);
        void WriteDensitiesPDB(string pdbfile);
        void Evaluate();
        void EvaluateFields();
        void EvaluateInduction();
        void EvaluateEnergy();
        
        bool Converged() { return _converged_R && _converged_K; }
        Property GenerateOutputString();
        string GenerateErrorString();
        
        virtual string IdentifyMethod() = 0;
        
        virtual EWD::triple<> ConvergeRealSpaceSum();
        virtual EWD::triple<> ConvergeReciprocalSpaceSum() = 0;
        virtual EWD::triple<> CalculateForegroundCorrection();
        virtual EWD::triple<> CalculateHigherRankCorrection();
        virtual EWD::triple<> CalculateShapeCorrection() { return EWD::triple<>(0.0,0.0,0.0); }
        virtual EWD::triple<> CalculateK0Correction() { return EWD::triple<>(0.0,0.0,0.0); }
        
        virtual void Field_ConvergeRealSpaceSum() { ; }
        virtual void Field_ConvergeReciprocalSpaceSum() { ; }
        virtual void Field_CalculateForegroundCorrection() { ; }
        virtual void Field_CalculateShapeCorrection() { ; }
        
        // To sort K-vectors via std::sort using a norm functor
        template<class Norm, class V>
        struct VectorSort
        {
            VectorSort() : _p(1e-40) { ; }
            VectorSort(double precision) : _p(precision) { ; }
            inline bool operator() (const V &v1, const V &v2);
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
            
        // K-vector class (for grading purposes)
        struct KVector
        {
            KVector(vec &k, double grade)
                : _k(k), _grade(grade) { ; }

            vec _k;
            double _grade;
            
            const vec &getK() const { return _k; }
            const double &getGrade() const { return _grade; }
            const double &getX() const { return _k.getX(); }
            const double &getY() const { return _k.getY(); }
            const double &getZ() const { return _k.getZ(); }            
        };
        
        // Specialized K-vector norm
        struct KNorm { inline double operator() (const KVector &v)
            { return -v.getGrade(); } };
        
    protected:
        
        EwdInteractor _ewdactor;
        XInteractor _actor;
        Logger *_log;
        
        // PERIODIC BOUNDARY
        Topology *_top;
        CSG::BoundaryCondition *_bc;    // Periodicity reduced to xy sub-space
        vec _center;
        
        // POLAR SEGMENTS
        // Part I - Ewald
        PolarTop *_ptop;
        vector< PolarSeg* > _bg_P;      // Period. density = _bg_N v _fg_N
        vector< PolarSeg* > _bg_N;      // Neutral background
        vector< PolarSeg* > _mg_N;      // Neutral midground
        vector< PolarSeg* > _fg_N;      // Neutral foreground
        vector< PolarSeg* > _fg_C;      // Charged foreground
        vector< bool > _inForeground;
        string _jobType;                // Calculated from _fg_C charge distr.
        // Part II - Thole
        vector< PolarSeg* > _polar_qm0;
        vector< PolarSeg* > _polar_mm1;
        vector< PolarSeg* > _polar_mm2; // Should not be used
        
        
        // CONVERGENCE
        // Part I - Ewald
        double _alpha;                  // _a = 1/(sqrt(2)*sigma)
        double _kfactor;
        double _rfactor;
        double _K_co;                   // k-space c/o
        double _R_co;                   // r-space c/o
        double _crit_dE;                // Energy convergence criterion [eV]
        bool   _converged_R;            // Did R-space sum converge?
        bool   _converged_K;            // Did K-space sum converge?
        bool   _field_converged_R;
        bool   _field_converged_K;
        // Part II - Thole
        bool _polar_do_induce;
        double _polar_aDamp;
        double _polar_wSOR_N;
        double _polar_wSOR_C;
        double _polar_cutoff;

        
        
        // LATTICE (REAL, RECIPROCAL)
        vec _a; vec _b; vec _c;         // Real-space lattice vectors
        int _na_max, _nb_max, _nc_max;  // Max. cell indices to sum over (R)
        vec _A; vec _B; vec _C;         // Reciprocal-space lattice vectors
        int _NA_max, _NB_max, _NC_max;  // Max. cell indices to sum over (K)
        double _LxLy;                   // |a^b|
        double _LxLyLz;                 // a*|b^c|
        
        VectorSort<MaxNorm,vec> _maxsort;
        VectorSort<EucNorm,vec> _eucsort;
        VectorSort<KNorm,KVector> _kvecsort;
        
        // ENERGIES
        // Part I - Ewald
        EWD::triple<> _ER;                     // R-space sum
        EWD::triple<> _EC;                     // R-space correction
        EWD::triple<> _EK;                     // K-space sum
        EWD::triple<> _E0;                     // K-space K=0 contribution
        EWD::triple<> _ET;                     // ER - EC + EK + E0
        EWD::triple<> _EDQ;                    // Higher-Rank FGC->MGN correction
        EWD::triple<> _EJ;                     // Geometry-dependent correction
        // Part II - Thole
        double _polar_ETT;
        double _polar_EPP;
        double _polar_EPU;
        double _polar_EUU;
        double _polar_EF00;
        double _polar_EF01;
        double _polar_EF02;
        double _polar_EF11;
        double _polar_EF12;
        double _polar_EM0;
        double _polar_EM1;
        double _polar_EM2;
        // I + II
        double _Estat;
        double _Eindu;
        double _Eppuu;
        
        
    };


template<class Norm, class V>
inline bool Ewald3DnD::VectorSort<Norm,V>::operator() (const V &v1,
    const V &v2) {
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