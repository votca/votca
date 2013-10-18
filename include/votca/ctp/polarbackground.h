#ifndef VOTCA_CTP_POLARBACKGROUND_H
#define VOTCA_CTP_POLARBACKGROUND_H

#include <votca/ctp/topology.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/ewdspace.h>
#include <votca/ctp/ewaldactor.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmthread.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace votca { namespace ctp {    
namespace EWD {

class PolarBackground
{
public:

    PolarBackground() : _top(NULL), _ptop(NULL), _log(NULL), _n_threads(1) {};
    PolarBackground(Topology *top, PolarTop *ptop, Property *opt, Logger *log, int n_threads);
   ~PolarBackground() {};
   
    void Threaded(int n_threads) { _n_threads = n_threads; }
    void Polarize();
    
    void FP_RealSpace();
    void FP_ReciprocalSpace();
    void FU_RealSpace(bool do_setup_nbs);
    void FU_ReciprocalSpace();
    
    // TODO Would be nicer to have as a local type in ::FP_RealSpace
    //      (which is possible by compiling with --std=c++0x)
    class FPThread : public QMThread
    {
    public:        
        FPThread(PolarBackground *master, int id);
       ~FPThread() { _full_bg_P.clear(); _part_bg_P.clear(); }       
        void Run(void);        
        void AddPolarSeg(PolarSeg *pseg) { _part_bg_P.push_back(pseg); return; }
        double Workload() { return 100.*_part_bg_P.size()/_full_bg_P.size(); }
        const int NotConverged() { return _not_converged_count; }
        double AvgRco() { return _avg_R_co; }
        void DoSetupNbs(bool do_setup) { _do_setup_nbs = do_setup; }
    private:        
        int _id;
        bool _verbose;
        PolarBackground *_master;
        EwdInteractor _ewdactor;
        vector<PolarSeg*> _full_bg_P;
        vector<PolarSeg*> _part_bg_P;
        int _not_converged_count;
        double _avg_R_co;
        bool _do_setup_nbs;
    };    
    
    class FUThread : public QMThread
    {
    public:        
        FUThread(PolarBackground *master, int id);
       ~FUThread() { _full_bg_P.clear(); _part_bg_P.clear(); }       
        void Run(void);        
        void AddPolarSeg(PolarSeg *pseg) { _part_bg_P.push_back(pseg); return; }        
        double Workload() { return 100.*_part_bg_P.size()/_full_bg_P.size(); }
        const int NotConverged() { return _not_converged_count; }
        double AvgRco() { return _avg_R_co; }
        void DoSetupNbs(bool do_setup) { _do_setup_nbs = do_setup; }
    private:
        int _id;
        bool _verbose;
        PolarBackground *_master;
        EwdInteractor _ewdactor;
        vector<PolarSeg*> _full_bg_P;
        vector<PolarSeg*> _part_bg_P;
        int _not_converged_count;
        double _avg_R_co;
        bool _do_setup_nbs;
    };
    
    
private:

    EwdInteractor _ewdactor;
    Logger *_log;
    int _n_threads;

    // PERIODIC BOUNDARY
    Topology *_top;
    string _shape;

    // POLAR SEGMENTS
    // Part I - Ewald
    PolarTop *_ptop;
    vector< PolarSeg* > _bg_P;      // Period. density = _bg_N v _fg_N

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
    double _polar_aDamp;
    double _polar_wSOR_N;
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

};


    
    
}
}}

#endif