#ifndef VOTCA_CTP_POLARBACKGROUND_H
#define VOTCA_CTP_POLARBACKGROUND_H

#include <votca/ctp/topology.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/ewdspace.h>
#include <votca/ctp/ewaldactor.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmthread.h>
#include <votca/ctp/threadforce.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace votca { namespace ctp {    
namespace EWD {

class PolarBackground
{
public:

    PolarBackground() : _top(NULL), _ptop(NULL), _log(NULL), _n_threads(1) {};
    PolarBackground(Topology *top, PolarTop *ptop, Property *opt, Logger *log, int n_threads);
   ~PolarBackground();
   
    void Threaded(int n_threads) { _n_threads = n_threads; }
    void Polarize();
    
    void FX_RealSpace(string mode, bool do_setup_nbs);    
    void FU_RealSpace(bool do_setup_nbs);
    void FX_ReciprocalSpace(string mode1, string mode2, bool generate_kvecs);
    
    void GenerateKVectors(vector<PolarSeg*> &ps1, vector<PolarSeg*> &ps2);
    
    // TODO Would be nicer to have RThread / KThread as a local type in ...
    //      ... but you cannot use local types as template parameters
    //      (which is however possible by compiling with --std=c++0x)
    
    class RThread : public QMThread
    {
    public:
        RThread(PolarBackground *master, bool do_setup_nbs) {
            _master = master;
            _do_setup_nbs = do_setup_nbs;
            _not_converged_count = 0;
            _ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);

            _mode_startfunct["FP_MODE"] = &RThread::FP_FieldCalc;
            _mode_startfunct["FU_MODE"] = &RThread::FU_FieldCalc;
            
            _mode_resetfunct["FP_MODE"] = &RThread::FX_FieldReset;
            _mode_resetfunct["FU_MODE"] = &RThread::FX_FieldReset;
            
            _mode_wloadfunct["FP_MODE"] = &RThread::FX_FieldWload;
            _mode_wloadfunct["FU_MODE"] = &RThread::FX_FieldWload;
        }
       ~RThread() {
           _full_bg_P.clear(); _part_bg_P.clear();
        }
       
        // Standard multi-mode thread interface
        typedef void (RThread::*StartFunct)();
        typedef void (RThread::*ResetFunct)();
        typedef double (RThread::*WloadFunct)();
        RThread *Clone() { return new RThread(_master, _do_setup_nbs); }
        const string &getMode() { return _current_mode; }
        void setMode(const string &mode) { _current_mode = mode; }
        void setVerbose(bool verbose) { _verbose = verbose; }        
        void Run(void) { 
            StartFunct start = _mode_startfunct[_current_mode]; 
            ((*this).*start)(); 
        }
        void Reset(string mode) {
            ResetFunct reset = _mode_resetfunct[mode];
            ((*this).*reset)();
        }
        double Workload(string mode) {
            WloadFunct wload = _mode_wloadfunct[mode];
            return ((*this).*wload)();
        }
        
        // Input
        void AddSharedInput(vector<PolarSeg*> &full_bg_P) { _full_bg_P = full_bg_P; }
        void AddAtomicInput(PolarSeg *pseg) { _part_bg_P.push_back(pseg); }
        
        // Mode targets
        void FU_FieldCalc();
        void FP_FieldCalc();
        void FX_FieldReset() { _not_converged_count = 0; _part_bg_P.clear(); }
        double FX_FieldWload() { return 100.*_part_bg_P.size()/_full_bg_P.size(); }
        
        // Convergence & output related
        const int NotConverged() { return _not_converged_count; }
        double AvgRco() { return _avg_R_co; }
        
    private:
        // Multi-mode thread members
        bool _verbose;
        string _current_mode;
        map<string,StartFunct> _mode_startfunct;
        map<string,ResetFunct> _mode_resetfunct;
        map<string,WloadFunct> _mode_wloadfunct;
        
        // Shared thread data
        PolarBackground *_master;
        vector<PolarSeg*> _full_bg_P;        
        // Atomized thread data
        EwdInteractor _ewdactor;
        vector<PolarSeg*> _part_bg_P;        
        // Convergence & output-related
        int _not_converged_count;
        double _avg_R_co;
        bool _do_setup_nbs;
    };    
    
    class KThread : public QMThread
    {
    public:
        
        KThread(PolarBackground *master) {
            _master = master;
            _ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);
            
            _mode_startfunct["SP_MODE"] = &KThread::SP_SFactorCalc;
            _mode_startfunct["FP_MODE"] = &KThread::FP_KFieldCalc;
            _mode_startfunct["SU_MODE"] = &KThread::SU_SFactorCalc;
            _mode_startfunct["FU_MODE"] = &KThread::FU_KFieldCalc;
            
            _mode_resetfunct["SP_MODE"] = &KThread::SFactorReset;
            _mode_resetfunct["FP_MODE"] = &KThread::KFieldReset;
            _mode_resetfunct["SU_MODE"] = &KThread::SFactorReset;
            _mode_resetfunct["FU_MODE"] = &KThread::KFieldReset;
            
            _mode_wloadfunct["SP_MODE"] = &KThread::SFactorWload;
            _mode_wloadfunct["FP_MODE"] = &KThread::KFieldWload;
            _mode_wloadfunct["SU_MODE"] = &KThread::SFactorWload;
            _mode_wloadfunct["FU_MODE"] = &KThread::KFieldWload;
        }
       ~KThread() {
           _full_bg_P.clear(); _full_kvecs.clear(); 
           _part_bg_P.clear(); _part_kvecs.clear();
        }
        
        // Standard multi-mode thread interface
        typedef void (KThread::*StartFunct)();
        typedef void (KThread::*ResetFunct)();
        typedef double (KThread::*WloadFunct)();
        KThread *Clone() { return new KThread(_master); }
        const string &getMode() { return _current_mode; }
        void setMode(const string &mode) { _current_mode = mode; }
        void setVerbose(bool verbose) { _verbose = verbose; }        
        void Run(void) { 
            StartFunct start = _mode_startfunct[_current_mode]; 
            ((*this).*start)(); 
        }
        void Reset(string mode) {
            ResetFunct reset = _mode_resetfunct[mode];
            ((*this).*reset)();
        }
        double Workload(string mode) {
            WloadFunct wload = _mode_wloadfunct[mode];
            return ((*this).*wload)();
        }    
        
        void AddSharedInput(vector<EWD::KVector*> &full_kvecs) { _full_kvecs.clear(); _full_kvecs = full_kvecs; }
        void AddSharedInput(vector<PolarSeg*> &full_bg_P) { _full_bg_P.clear(); _full_bg_P = full_bg_P; }
        
        // MODE 1 : Compute structure factor of _part_kvecs using _full_bg_P
        void SP_SFactorCalc();
        void SU_SFactorCalc();
        void AddAtomicInput(EWD::KVector *k) { _part_kvecs.push_back(k); }
        void SFactorReset() { _part_kvecs.clear(); _full_kvecs.clear(); }
        double SFactorWload() { return (double)_part_kvecs.size()/_full_kvecs.size(); }
        
        // MODE 2 : Increment fields of _part_bg_P using _full_kvecs
        void FP_KFieldCalc();
        void FU_KFieldCalc();
        void AddAtomicInput(PolarSeg *pseg) { _part_bg_P.push_back(pseg); }
        void KFieldReset() { ; }
        double KFieldWload() { return (double)_part_bg_P.size()/_full_bg_P.size(); }
        
    private:
        // Multi-mode thread members
        bool _verbose;
        string _current_mode;
        map<string,StartFunct> _mode_startfunct;
        map<string,ResetFunct> _mode_resetfunct;
        map<string,WloadFunct> _mode_wloadfunct;
        
        // Shared thread data
        PolarBackground *_master;
        vector<PolarSeg*> _full_bg_P;
        vector<EWD::KVector*> _full_kvecs;
        
        // Atomized thread data
        EwdInteractor _ewdactor;
        vector<EWD::KVector*> _part_kvecs;
        vector<PolarSeg*> _part_bg_P;
    
    public:
        // Convergence info
        double _rms_sum_re;
        double _sum_im;
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

    vector<EWD::KVector*> _kvecs_2_0;   // K-vectors with two components = zero
    vector<EWD::KVector*> _kvecs_1_0;   // K-vectors with one component  = zero
    vector<EWD::KVector*> _kvecs_0_0;   // K-vectors with no  components = zero
    double _kxyz_s1s2_norm;

};    

}}}

#endif