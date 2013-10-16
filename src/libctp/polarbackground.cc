#include <votca/ctp/polarbackground.h>
#include <boost/format.hpp>

namespace votca {
namespace ctp {
namespace EWD {

using boost::format;
    
PolarBackground::PolarBackground(Topology *top, PolarTop *ptop, Property *opt, 
    Logger *log) : _top(top), _ptop(ptop), _log(log) {
    
    // EVALUATE OPTIONS
    string pfx = "options.ewdbgpol";
    // Ewald parameters
    _R_co = opt->get(pfx+".coulombmethod.cutoff").as<double>();
    _crit_dE = opt->get(pfx+".convergence.energy").as<double>();
    if (opt->exists(pfx+".convergence.kfactor"))
        _kfactor = opt->get(pfx+".convergence.kfactor").as<double>();
    else
        _kfactor = 100.;
    if (opt->exists(pfx+".convergence.rfactor"))
        _rfactor = opt->get(pfx+".convergence.rfactor").as<double>();
    else
        _rfactor = 6.;    
    // Polar parameters    
    if (opt->exists(pfx+".polarmethod.cutoff")) 
        _polar_cutoff = opt->get(pfx+".polarmethod.cutoff").as<double>();
    else
        _polar_cutoff = 0.0;
    if (opt->exists(pfx+".polarmethod.wSOR_N"))
        _polar_wSOR_N = opt->get(pfx+".polarmethod.wSOR_N").as<double>();
    else
        _polar_wSOR_N = 0.35;
    if (opt->exists(pfx+".polarmethod.aDamp"))
        _polar_aDamp = opt->get(pfx+".polarmethod.aDamp").as<double>();
    else
        _polar_aDamp = 0.390;
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _K_co = _kfactor/_R_co;
    _alpha = _rfactor/_R_co;
    _ewdactor = EwdInteractor(_alpha, _polar_aDamp);
    
    // SET-UP REAL & RECIPROCAL SPACE
    _a = _top->getBox().getCol(0);
    _b = _top->getBox().getCol(1);
    _c = _top->getBox().getCol(2);
    _LxLyLz = _a*(_b^_c);
    _LxLy = abs(_a ^ _b);
    
    _A = 2*M_PI/_LxLyLz * _b^_c;
    _B = 2*M_PI/_LxLyLz * _c^_a;
    _C = 2*M_PI/_LxLyLz * _a^_b;

    _na_max = ceil(_R_co/maxnorm(_a)-0.5)+1;
    _nb_max = ceil(_R_co/maxnorm(_b)-0.5)+1;
    _nc_max = ceil(_R_co/maxnorm(_c)-0.5)+1;

    _NA_max = ceil(_K_co/maxnorm(_A));
    _NB_max = ceil(_K_co/maxnorm(_B));
    _NC_max = ceil(_K_co/maxnorm(_C));
    
    // SET-UP POLAR GROUNDS (FORE-, MID-, BACK-)
    _bg_P.clear();
    _bg_P = ptop->BGN();
    
    // CALCULATE COG POSITIONS, NET CHARGE
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    double Q_bg_P = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }
    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size:")).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush;
    
    if (std::abs(Q_bg_P) > 1e-2) {
        cout << endl;
        cout << endl << format("***************************** ERROR ******************************");
        cout << endl << format("       Background charge |Q(BGP)| is larger than 0.01e.");
        cout << endl << format("       Be more precise: e.g. rounding error?");
        cout << endl << format("       Or think again: e.g. erroneous parametrization?");
        cout << endl << format("******************************************************************");
        cout << endl;
    }
    
    // CHARGE APPROPRIATELY & DEPOLARIZE
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    
    // CALCULATE NET DIPOLE OF BGP & FGC
    vec netdpl_bgP = vec(0,0,0);
    double qzz_bgP = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_bgP += (*pit)->getPos() * (*pit)->getQ00();
            qzz_bgP += (*pit)->getQ00() * ((*pit)->getPos().getZ() * (*pit)->getPos().getZ());
        }
    }
    
    LOG(logINFO,*_log)
        << (format("Net dipole moment of background density")).str()
        << flush << (format("  o D(BGP) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_bgP.getX() % netdpl_bgP.getY() % netdpl_bgP.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o Sigma q|z|**2 [e*nm**2] = %1$+1.7f   ")
        % qzz_bgP) << flush;
    
    return;
}
    
    
    
    
    
    
    
}}}