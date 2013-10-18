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
    void FU_RealSpace();
    void FU_ReciprocalSpace();
    
    // TODO Would be nicer to have as a local type in ::FP_RealSpace
    //      ... but this is only possible via --std=c++0x
    class FPThread : public QMThread
    {
    public:
        
        FPThread(PolarBackground *master, int id) : _id(id) {
            _master = master;
            _full_bg_P = master->_bg_P;
            _ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);
            _verbose = (_id-1 == (_full_bg_P.size()-1) % _master->_n_threads);
        }
       ~FPThread() {
            _full_bg_P.clear(); 
            _part_bg_P.clear(); 
        }
       
        void Run(void) {
            vector<PolarSeg*>::iterator sit1; 
            vector<APolarSite*> ::iterator pit1;
            vector<PolarSeg*>::iterator sit2; 
            vector<APolarSite*> ::iterator pit2;
            double rms = 0.0;
            int rms_count = 0;
            for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
                PolarSeg *pseg1 = *sit1;
                if (_verbose)
                    LOG(logDEBUG,*(_master->_log))
                        << "\rMST DBG     - Progress " << pseg1->getId() 
                        << "/" << _full_bg_P.size() << flush;
                for (sit2 = _full_bg_P.begin(); sit2 < _full_bg_P.end(); ++sit2) {
                    PolarSeg *pseg2 = *sit2;
                    // Identical?
                    if (pseg1 == pseg2) continue;
                    if (pseg1->getId() == pseg2->getId()) assert(false);
                    // Apply periodic-boundary correction, check c/o, shift
                    vec dr12_pbc = _master->_top->PbShortestConnect(pseg1->getPos(), pseg2->getPos());
                    if (votca::tools::abs(dr12_pbc) > _master->_R_co) continue;
                    vec dr12_dir = pseg2->getPos() - pseg1->getPos();
                    vec s22x = dr12_pbc - dr12_dir;
                    // Interact taking into account shift
                    for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
                        for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
                            rms += _ewdactor.FP12_ERFC_At_By(*(*pit1), *(*pit2), s22x);
                            rms_count += 1;
                        }
                    }
                }
            }
            rms = sqrt(rms/rms_count)*EWD::int2V_m;
            return;
        }
        
        void AddPolarSeg(PolarSeg *pseg) {
            _part_bg_P.push_back(pseg);
            return;
        }
        
        double Workload() { return 100.*(double)_part_bg_P.size()/_full_bg_P.size(); }
        
    private:
        
        int _id;
        bool _verbose;
        PolarBackground *_master;
        EwdInteractor _ewdactor;
        vector<PolarSeg*> _full_bg_P;
        vector<PolarSeg*> _part_bg_P;
    };
    
    
    class FUThread : public QMThread
    {
    public:
        
        FUThread(PolarBackground *master, int id) : _id(id) {
            _master = master;
            _full_bg_P = master->_bg_P;
            _ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);
            _verbose = (_id-1 == (_full_bg_P.size()-1) % _master->_n_threads);
        }
       ~FUThread() {
            _full_bg_P.clear(); 
            _part_bg_P.clear(); 
        }
       
        void Run(void) {
            vector<PolarSeg*>::iterator sit1; 
            vector<APolarSite*> ::iterator pit1;
            vector<PolarSeg*>::iterator sit2; 
            vector<APolarSite*> ::iterator pit2;
            double rms = 0.0;
            int rms_count = 0;
            for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
                PolarSeg *pseg1 = *sit1;
                if (_verbose)
                    LOG(logDEBUG,*(_master->_log))
                        << "\rMST DBG     - Progress " << pseg1->getId() 
                        << "/" << _full_bg_P.size() << flush;
                for (sit2 = _full_bg_P.begin(); sit2 < _full_bg_P.end(); ++sit2) {
                    PolarSeg *pseg2 = *sit2;
                    // Identical?
                    if (pseg1 == pseg2) continue;
                    if (pseg1->getId() == pseg2->getId()) assert(false);
                    // Apply periodic-boundary correction, check c/o, shift
                    vec dr12_pbc = _master->_top->PbShortestConnect(pseg1->getPos(), pseg2->getPos());
                    if (votca::tools::abs(dr12_pbc) > _master->_R_co) continue;
                    vec dr12_dir = pseg2->getPos() - pseg1->getPos();
                    vec s22x = dr12_pbc - dr12_dir;
                    // Interact taking into account shift
                    for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
                        for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
                            rms += _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2), s22x);
                            rms_count += 1;
                        }
                    }
                }
            }
            rms = sqrt(rms/rms_count)*EWD::int2V_m;
            return;
        }
        
        void AddPolarSeg(PolarSeg *pseg) {
            _part_bg_P.push_back(pseg);
            return;
        }
        
        double Workload() { return 100.*(double)_part_bg_P.size()/_full_bg_P.size(); }
        
    private:
        
        int _id;
        bool _verbose;
        PolarBackground *_master;
        EwdInteractor _ewdactor;
        vector<PolarSeg*> _full_bg_P;
        vector<PolarSeg*> _part_bg_P;
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