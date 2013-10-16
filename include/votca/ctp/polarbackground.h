#ifndef VOTCA_CTP_POLARBACKGROUND_H
#define VOTCA_CTP_POLARBACKGROUND_H

#include <votca/ctp/topology.h>
#include <votca/ctp/polartop.h>
#include <votca/ctp/ewdspace.h>
#include <votca/ctp/ewaldactor.h>
#include <votca/ctp/logger.h>

namespace votca { namespace ctp {    
namespace EWD {

class PolarBackground
{
public:

    PolarBackground(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
   ~PolarBackground() {};
    
private:

    EwdInteractor _ewdactor;
    Logger *_log;

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