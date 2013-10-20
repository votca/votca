#include <votca/ctp/polarbackground.h>
#include <boost/format.hpp>

namespace votca {
namespace ctp {
namespace EWD {

using boost::format;
    
PolarBackground::PolarBackground(Topology *top, PolarTop *ptop, Property *opt, 
    Logger *log, int n_threads) : _top(top), _ptop(ptop), _log(log), 
    _n_threads(n_threads) {
    
    // EVALUATE OPTIONS
    string pfx = "options.ewdbgpol";
    // Ewald parameters
    _shape = opt->get(pfx+".coulombmethod.shape").as<string>();
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


void PolarBackground::Polarize() {    
    
    TLogLevel dbg = logDEBUG;
    TLogLevel inf = logINFO;
    TLogLevel err = logERROR;
    Logger &log = *_log;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    /*
    Verify neutrality & depolarize
    Generate permanent fields (FP)
      o Converge intermolecular real-space contribution, remember cut-off
      o Converge reciprocal-space contribution, remember K-vectors
      o Calculate shape fields
      o Apply MOLECULAR foreground correction     
    Induce to 1st order
    Loop until 2nd-order fields converged
      | Reset 2nd-order fields
      | (Re-)generate induction fields (FU)
      | o Real-space INTRAmolecular contribution to 2nd-order fields
      | o Real-space INTERmolecular contribution to 2nd-order fields
      | o Reciprocal-space contribution, work off remembered K-vectors
      | o Calculate shape fields
      | o Apply ATOMIC foreground correction
      | Induce to 2nd order
      + Check convergence
    Extract (or serialize) induction state to hard-drive
    */
    
    // I GENERATE PERMANENT FIELDS (FP)
    LOG(dbg,log) << flush;
    LOG(dbg,log) << "Generate permanent fields (FP)" << flush;
    // I.A Intermolecular real-space contribution
    this->FP_RealSpace();
    // I.B Reciprocal-space contribution
    LOG(dbg,log) << "  o TODO Reciprocal-space" << flush;
    ;
    // I.C Shape fields
    LOG(dbg,log) << "  o Shape fields" << flush;
    if (_shape == "xyslab") {
        double TwoPi_V = 2*M_PI/_LxLyLz;        
        _ewdactor.FP12_XYSlab_ShapeField_At_By(_bg_P, _bg_P, TwoPi_V);        
    }
    else {
        LOG(logERROR,*_log)
            << (format("Shape '%1$s' not implemented. Omitting shape fields.") 
            % _shape) << flush;
    }
    // I.D Molecular ERF self-interaction correction
    LOG(dbg,log) << "  o Molecular SI correction" << flush;
    double rms = 0.0;
    int rms_count = 0;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            for (pit2 = (*sit1)->begin(); pit2 < (*sit1)->end(); ++pit2) {
                rms += _ewdactor.FP12_ERF_At_By(*(*pit1), *(*pit2));
                rms_count += 1;
            }
        }
    }
    rms = sqrt(rms/rms_count)*EWD::int2V_m;
    
    
    // II INDUCE TO 1ST ORDER
    LOG(dbg,log) << "Induce to first order" << flush;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            (*pit1)->InduceDirect();
        }
    }
    
    // III CONVERGE INDUCTION FIELDS (FU)
    int iter = 0;
    int max_iter = 512;
    double epstol = 1e-3;
    for ( ; iter < max_iter; ++iter) {
        LOG(dbg,log) << flush;
        LOG(dbg,log) << "Iter " << iter << " started" << flush;
        
        // III.A Reset 2nd-order fields (FU)
        LOG(dbg,log) << "  o Reset fields (FU)" << flush;
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                (*pit1)->ResetFieldU();
            }
        }
        // III.B (Re-)generate induction fields FU
        LOG(dbg,log) << "  o (Re-)generate induction fields" << flush;
        // (1) Real-space intramolecular contribution
        LOG(dbg,log) << "  o Real-space, intramolecular" << flush;
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = pit1+1; pit2 < (*sit1)->end(); ++pit2) {
                    _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2));
                    _ewdactor.FU12_ERFC_At_By(*(*pit2), *(*pit1));
                }
            }
        }
        // (2) Real-space intermolecular contribution
        this->FU_RealSpace();
        // (3) Reciprocal-space contribution
        LOG(dbg,log) << "  o TODO Reciprocal-space" << flush;
        ;
        // (4) Calculate shape fields
        LOG(dbg,log) << "  o Shape fields" << flush;
        if (_shape == "xyslab") {
            double TwoPi_V = 2*M_PI/_LxLyLz;        
            _ewdactor.FU12_XYSlab_ShapeField_At_By(_bg_P, _bg_P, TwoPi_V);        
        }
        else {
            LOG(logERROR,*_log) 
                << (format("Shape '%1$s' not implemented. Omitting shape fields.") 
                % _shape) << flush;
        }
        // (5) Apply atomic ERF self-interaction correction
        LOG(dbg,log) << "    - Atomic SI correction" << flush;
        rms = 0.0;
        rms_count = 0;
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                rms += _ewdactor.FU12_ERF_At_By(*(*pit1), *(*pit1));
                rms_count += 1;
            }
        }
        rms = sqrt(rms/rms_count)*EWD::int2V_m;
        
        // III.C Induce again
        LOG(dbg,log) << "  o Induce again" << flush;
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                (*pit1)->Induce(_polar_wSOR_N);
            }
        }
        
        // III.D Check for convergence
        LOG(dbg,log) << "  o Convergence check" << flush;
        bool converged = true;
        double maxdU = -1;
        double avgdU = 0.0;
        int baseN = 0;
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                double dU = (*pit1)->HistdU();
                avgdU += dU;
                ++baseN;
                if (dU > maxdU) { maxdU = dU; }
                if (dU > epstol) { converged = false; }
            }
        }
        avgdU /= baseN;
        if (avgdU < epstol*0.1) { converged = true; }
        if (converged) {
            LOG(dbg,log) << flush;
            LOG(dbg,log) << ":: Converged induction fields" << flush;
            break;
        }
        else if (iter == max_iter-1) {
            throw std::runtime_error("Not converged.");
            break;
        }
    }
        
    return;
}


void PolarBackground::FP_RealSpace() {
    
    TLogLevel dbg = logDEBUG;
    TLogLevel inf = logINFO;
    TLogLevel err = logERROR;
    Logger &log = *_log;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    LOG(dbg,log) << 
        (format("  o Real-space intermolecular FP R(co)=%1$1.1fnm, N(th)=%2$d") 
        % _R_co % _n_threads) << flush;
    
    // Create threads
    vector<FPThread*> fpthreads;
    for (int t = 0; t < _n_threads; ++t) {
        FPThread *newthread = new FPThread(this, t+1);
        fpthreads.push_back(newthread);
    }
    
    // Distribute workload
    LOG(dbg,log) << "    - Thread workload = [ ";
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        int tidx = ((*sit1)->getId()-1) % _n_threads;
        fpthreads[tidx]->AddPolarSeg(*sit1);
    }
    for (int t = 0; t < _n_threads; ++t) {
        LOG(dbg,log) << (format("%1$1.2f%% ") % fpthreads[t]->Workload());
    }
    LOG(dbg,log) << "]" << flush;
    
    // Start & wait
    LOG(dbg,log) << "    - Start & wait until done" << flush << flush;
    _log->setPreface(logDEBUG, "");
    for (int t = 0; t < _n_threads; ++t) fpthreads[t]->Start();
    for (int t = 0; t < _n_threads; ++t) fpthreads[t]->WaitDone();
    for (int t = 0; t < _n_threads; ++t) delete fpthreads[t];
    _log->setPreface(logDEBUG,   "\nMST DBG");
    
    
//    cout << endl;
//    double rms = 0.0;
//    int rms_count = 0;
//    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
//        PolarSeg *pseg1 = *sit1;
//        cout << "\r " << pseg1->getId() << flush;
//        for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
//            PolarSeg *pseg2 = *sit2;
//            // Identical?
//            if (pseg1 == pseg2) continue;
//            if (pseg1->getId() == pseg2->getId()) assert(false);
//            // Apply periodic-boundary correction, check c/o, shift
//            vec dr12_pbc = _top->PbShortestConnect(pseg1->getPos(), pseg2->getPos());
//            if (votca::tools::abs(dr12_pbc) > _R_co) continue;
//            vec dr12_dir = pseg2->getPos() - pseg1->getPos();
//            vec s22x = dr12_pbc - dr12_dir;
//            // Interact taking into account shift
//            for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
//                for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
//                    rms += _ewdactor.FP12_ERFC_At_By(*(*pit1), *(*pit2), s22x);
//                    rms_count += 1;
//                }
//            }
//        }
//    }
//    rms = sqrt(rms/rms_count)*EWD::int2V_m;
    
    return;
}


void PolarBackground::FU_RealSpace() {
    
    TLogLevel dbg = logDEBUG;
    TLogLevel inf = logINFO;
    TLogLevel err = logERROR;
    Logger &log = *_log;
    
    LOG(dbg,log) << 
        (format("  o Real-space intermolecular FU R(co)=%1$1.1fnm, N(th)=%2$d") 
        % _R_co % _n_threads) << flush;
    
    // Create threads
    vector<FUThread*> futhreads;
    for (int t = 0; t < _n_threads; ++t) {
        FUThread *newthread = new FUThread(this, t+1);
        futhreads.push_back(newthread);
    }
    
    // Distribute workload
    LOG(dbg,log) << "    - Thread workload = [ ";
    vector<PolarSeg*>::iterator sit1;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        int tidx = ((*sit1)->getId()-1) % _n_threads;
        futhreads[tidx]->AddPolarSeg(*sit1);
    }
    for (int t = 0; t < _n_threads; ++t) {
        LOG(dbg,log) << (format("%1$1.2f%% ") % futhreads[t]->Workload());
    }
    LOG(dbg,log) << "]" << flush;
    
    // Start & wait
    LOG(dbg,log) << "    - Start & wait until done" << flush << flush;
    _log->setPreface(logDEBUG, "");
    for (int t = 0; t < _n_threads; ++t) futhreads[t]->Start();
    for (int t = 0; t < _n_threads; ++t) futhreads[t]->WaitDone();
    for (int t = 0; t < _n_threads; ++t) delete futhreads[t];
    _log->setPreface(logDEBUG,   "\nMST DBG");
    
    return;
}

    
}}}