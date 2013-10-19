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
        virtual void GenerateKVectors(vector<PolarSeg*> &ps1, vector<PolarSeg*> &ps2) { ; }
        
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
        
        virtual void PolarizeBackground() { ; }
        
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
        
        
        // TASKS
        bool _task_polarize_bg;
        bool _task_calculate_fields;
        bool _task_polarize_fg;
        bool _task_evaluate_energy;
        
        
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
        bool   _did_field_pin_R_shell_idx;
        int    _field_R_shell_idx;
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
        
        EWD::VectorSort<EWD::MaxNorm,vec> _maxsort;
        EWD::VectorSort<EWD::EucNorm,vec> _eucsort;
        EWD::VectorSort<EWD::KNorm,EWD::KVector> _kvecsort;
        
        vector<EWD::KVector> _kvecs_2_0;     // K-vectors with two components = zero
        vector<EWD::KVector> _kvecs_1_0;     // K-vectors with one component  = zero
        vector<EWD::KVector> _kvecs_0_0;     // K-vectors with no  components = zero
        double _kxyz_s1s2_norm;
        bool _did_generate_kvectors;
        
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


}}


#endif