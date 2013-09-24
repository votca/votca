#include <votca/ctp/ewaldnd.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/math/special_functions/round.hpp>


using boost::format;


namespace votca { namespace ctp {


Ewald3DnD::~Ewald3DnD() {
    vector< PolarSeg* >::iterator sit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _polar_qm0.clear();
    _polar_mm1.clear();
    _polar_mm2.clear();
}
    
    
Ewald3DnD::Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
    : _top(top), _ptop(ptop), _log(log) {
    
    // EVALUATE OPTIONS
    string pfx = "options.ewald";
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
        _rfactor = 3.5;    
    // Polar parameters    
    if (opt->exists(pfx+".polarmethod.induce"))
        _polar_do_induce = opt->get(pfx+".polarmethod.induce").as<bool>();
    else
        _polar_do_induce = false;
    if (opt->exists(pfx+".polarmethod.cutoff")) 
        _polar_cutoff = opt->get(pfx+".polarmethod.cutoff").as<double>();
    else
        _polar_cutoff = 0.0;
    if (opt->exists(pfx+".polarmethod.wSOR_N"))
        _polar_wSOR_N = opt->get(pfx+".polarmethod.wSOR_N").as<double>();
    else
        _polar_wSOR_N = 0.35;
    if (opt->exists(pfx+".polarmethod.wSOR_C"))
        _polar_wSOR_C = opt->get(pfx+".polarmethod.wSOR_C").as<double>();
    else
        _polar_wSOR_C = 0.30;
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
    _center = ptop->getCenter();
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _fg_C = ptop->FGC();
    _fg_N = ptop->FGN();
    _bg_N = ptop->BGN();
    
    // Grow fg vs. bg according to induction cut-off
    this->ExpandForegroundReduceBackground(_polar_cutoff);    
    _bg_P.insert(_bg_P.end(), _fg_N.begin(), _fg_N.end());
    _bg_P.insert(_bg_P.end(), _bg_N.begin(), _bg_N.end());    
    _inForeground.resize(_bg_P.size()+1,false);
    
    // SET-UP MIDGROUND (INCLUDING PERIODIC IMAGES IF REQUIRED)
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Generate periodic images. ";
    this->SetupMidground(_R_co);
    
    // CALCULATE COG POSITIONS, NET CHARGE; SET-UP BOOLEAN FG TABLE
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    double Q_fg_C = 0.0;
    double Q_fg_C_2nd = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        double Qseg = (*sit)->CalcTotQ();
        Q_fg_C += Qseg;
        Q_fg_C_2nd += Qseg*Qseg / _fg_C.size();
        _inForeground[(*sit)->getId()] = true;
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_N += (*sit)->CalcTotQ();
        assert(_inForeground[(*sit)->getId()] == true);
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_N += (*sit)->CalcTotQ();
        assert(_inForeground[(*sit)->getId()] == false);
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }
    
    // DETERMINE JOB TYPE
    int iQ1 = boost::math::iround(Q_fg_C);
    int iQ2 = boost::math::iround(Q_fg_C_2nd);
    if (iQ1 == +1)         _jobType = "hole-like";
    else if (iQ1 == -1)    _jobType = "electron-like";
    else if (iQ1 == 0 && iQ2 == 0) _jobType = "neutral";
    else if (iQ1 == 0 && iQ2 > 0) _jobType = "charge-transfer-like";
    else _jobType = "bipolaron-like";

    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size:")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| ~ %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush << (format("  o Job type '%3$s' (iQ1=%1$d, iQ2=%2$d)") % iQ1 % iQ2 % _jobType).str()
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
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    
    // CALCULATE NET DIPOLE OF BGP & FGC
    vec netdpl_bgP = vec(0,0,0);
    vec netdpl_fgC = vec(0,0,0);
    double qzz_bgP = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_bgP += (*pit)->getPos() * (*pit)->getQ00();
            qzz_bgP += (*pit)->getQ00() * ((*pit)->getPos().getZ() * (*pit)->getPos().getZ());
        }
    }
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_fgC += (*pit)->getPos() * (*pit)->getQ00();
        }
    }
    
    LOG(logINFO,*_log)
        << (format("Net dipole moment of background density")).str()
        << flush << (format("  o D(BGP) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_bgP.getX() % netdpl_bgP.getY() % netdpl_bgP.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o D(FGC) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_fgC.getX() % netdpl_fgC.getY() % netdpl_fgC.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o Sigma q|z|**2 [e*nm**2] = %1$+1.7f   ")
        % qzz_bgP) << flush;
    
    return;
}


void Ewald3DnD::ExpandForegroundReduceBackground(double polar_R_co) {
    
    vector<PolarSeg*>::iterator sit1;
    vector<PolarSeg*>::iterator sit2;
    
    for (sit1 = _polar_mm1.begin(); sit1 < _polar_mm1.end(); ++sit1)
        delete *sit1;
    _polar_qm0.clear();
    _polar_mm1.clear();
    _polar_mm2.clear();
    
    // Target containers
    vector<PolarSeg*> new_exp_fg_C;
    vector<PolarSeg*> exp_fg_N;
    vector<PolarSeg*> red_bg_N;    
    
    // Foreground remains in foreground + contributes to QM0
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        new_exp_fg_C.push_back(*sit1);
        _polar_qm0.push_back(*sit1);
    }
    for (sit1 = _fg_N.begin(); sit1 < _fg_N.end(); ++sit1) {
        exp_fg_N.push_back(*sit1);
    }
    
    // Shift segments from back- to foreground based on distance vs. c/o
    // 'new' foreground copies contribute to MM1
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _bg_N.begin(); sit2 < _bg_N.end(); ++sit2) {
            PolarSeg *seg1 = *sit1;
            PolarSeg *seg2 = *sit2;            
            double dR = votca::tools::abs(seg1->getPos() - seg2->getPos());
            if (dR < polar_R_co) {
                // 'new' copy required for induction that affects fgC,
                // ... but NOT fgN
                PolarSeg *fgCopy = new PolarSeg(seg2);
                new_exp_fg_C.push_back(fgCopy);
                exp_fg_N.push_back(seg2);
                _polar_mm1.push_back(fgCopy);
            }
            else {
                red_bg_N.push_back(seg2);
            }
        }
    }
    
    // Exchange new for old containers
    _fg_C.clear();
    _fg_N.clear();
    _bg_N.clear();
    _fg_C = new_exp_fg_C;
    _fg_N = exp_fg_N;
    _bg_N = red_bg_N;
    
    bool clean = false; // Already deleted via fgC
    _ptop->setQM0(_polar_qm0, clean);
    _ptop->setMM1(_polar_mm1, clean);
    _ptop->setMM2(_polar_mm2, clean);
    
    assert(_polar_qm0.size()+_polar_mm1.size() == _fg_C.size());
    assert(_fg_C.size() == _fg_N.size());
    
    return;
}


void Ewald3DnD::SetupMidground(double R_co) {
    // SET-UP MIDGROUND
    // TODO Extend this to several molecules in the foreground - DONE?
    // NOTE No periodic-boundary correction here: We require that all
    //      polar-segment CoM coords. be folded with respect to the central
    //      segment
    // NOTE Excludes interaction of polar segment with neutral self in real-
    //      space sum
    // NOTE Includes periodic images if within cut-off
    
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _mg_N.clear();
    
    // IMAGE BOXES TO CONSIDER
    vector< int > nas;
    vector< int > nbs;
    vector< int > ncs;
    for (int na = -_na_max; na < _na_max+1; ++na)
        nas.push_back(na);
    for (int nb = -_nb_max; nb < _nb_max+1; ++nb)
        nbs.push_back(nb);
    for (int nc = -_nc_max; nc < _nc_max+1; ++nc)
        ncs.push_back(nc);
    vector< int >::iterator ait;
    vector< int >::iterator bit;
    vector< int >::iterator cit;
    
    // SAMPLE MIDGROUND FROM BGP EXCLUDING CENTRAL SEG.
    assert(_fg_N.size() == _fg_C.size());
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        Segment *other = _top->getSegment(pseg->getId());
        // Periodic images
        for (ait = nas.begin(); ait < nas.end(); ++ait) {
            int na = *ait;
            for (bit = nbs.begin(); bit < nbs.end(); ++bit) {                
                int nb = *bit;
                for (cit = ncs.begin(); cit < ncs.end(); ++cit) {
                    int nc = *cit;
                    if (na == 0 && nb == 0 && nc == 0 && _inForeground[other->getId()])
                        continue;
                    vec L = na*_a + nb*_b + nc*_c;
                    vec pos_L = pseg->getPos() + L;
                    double dR_L = abs(_center-pos_L);
                    if (dR_L <= R_co) {
                        PolarSeg *newSeg = new PolarSeg(pseg);
                        newSeg->Translate(L);
                        _mg_N.push_back(newSeg);
                    }
                } // Loop over nc
            } // Loop over nb
        } // Loop over na
    } // Loop over BGP
    return;
}


void Ewald3DnD::WriteDensitiesPDB(string pdbfile) {
    // COORDINATE OUTPUT FOR VISUAL CHECK
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;    
    FILE *out;
    out = fopen(pdbfile.c_str(),"w");
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGC");
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGN");
        }
    }
    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "QM0");
        }
    }
    for (sit = _polar_mm1.begin(); sit < _polar_mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM1");
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MGN");
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGN");
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGP");
        }
    }
    fclose(out);
    return;
}


void Ewald3DnD::Evaluate() {
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "System & Ewald parameters (" << IdentifyMethod() << ")" << flush;
    LOG(logDEBUG,*_log) << "  o Real-space unit cell:      " << _a << " x " << _b << " x " << _c << flush;
    LOG(logDEBUG,*_log) << "  o Real-space c/o (guess):    " << _R_co << " nm" << flush;
    LOG(logDEBUG,*_log) << "  o na(max), nb(max), nc(max): " << _na_max << ", " << _nb_max << ", " << _nc_max << flush;
    LOG(logDEBUG,*_log) << "  o 1st Brillouin zone:        " << _A << " x " << _B << " x " << _C << flush;
    LOG(logDEBUG,*_log) << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o R-K switching param.       " << _alpha << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*_log) << "  o LxLy (for 3D2D EW):        " << _LxLy << " nm**2" << flush;
    LOG(logDEBUG,*_log) << "  o kx(max), ky(max), kz(max): " << _NA_max << ", " << _NB_max << ", " << _NC_max << flush;
    
    EvaluateFields();
    EvaluateInduction();
    EvaluateEnergy();
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN)  = %1$+1.7f eV") % _ER).str()
        << flush << (format("  + EKK(FGC->BGP)  = %1$+1.7f eV") % _EK).str()       
        << flush << (format("  - EPP(FGC->FGN)  = %1$+1.7f eV") % _EC).str()
        << flush << (format("  = ------------------------------")).str()
        << flush << (format("  + SUM(E)         = %1$+1.7f eV") % (_ET-_EJ-_EDQ-_E0)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + EK0(FGC->BGP)  = %1$+1.7f eV") % _E0).str() 
        << flush << (format("  + EDQ(FGC->MGN)  = %1$+1.7f eV") % _EDQ).str()
        << flush << (format("  + EJ(shape-dep.) = %1$+1.7f eV") % _EJ).str()
        << flush << (format("  = ------------------------------")).str()
        << flush << (format("  + SUM(E) (0,Q,J) = %1$+1.7f eV") % (_ET)).str()
        << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC <> FGC")).str()
        << flush << (format("  + EPP(FGC<>FGC)  = %1$+1.7f eV") % _polar_EPP).str()
        << flush << (format("  + EPU(FGC<>FGC)  = %1$+1.7f eV") % _polar_EPU).str()       
        << flush << (format("  + EUU(FGC<>FGC)  = %1$+1.7f eV") % _polar_EUU).str()
        << flush << (format("  = ------------------------------")).str()
        << flush << (format("  + SUM(E)         = %1$+1.7f eV") % _polar_ETT).str()
        << flush << (format("    ==============================")).str()
        << flush << (format("    SUM(E) (1,2)   = %1$+1.7f eV") % (_ET+_polar_ETT)).str()
        << flush;
    LOG(logDEBUG,*_log) << flush;
    return;
}


void Ewald3DnD::EvaluateFields() {
    
    // RESET PERMANENT FIELDS
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }

    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeRealSpaceSum();    

    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeReciprocalSpaceSum();

    // SHAPE-CORRECTION (3D3D)/ K0-CORRECTION (3D2D)
    Field_CalculateShapeCorrection();

    // FOREGROUND CORRECTION (3D2D && 3D3D)
    Field_CalculateForegroundCorrection();
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    LOG(logDEBUG,*_log) << flush << "Foreground fields:" << flush;
    int fieldCount = 0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {        
        PolarSeg* pseg = *sit1;        
        for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            LOG(logDEBUG,*_log)
               << (format("F = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m") 
                    % (fp.getX()*_ewdactor.int2V_m) 
                    % (fp.getY()*_ewdactor.int2V_m) 
                    % (fp.getZ()*_ewdactor.int2V_m)).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*_log)
                    << "F = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
    
    return;
}


void Ewald3DnD::EvaluateInduction() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << format("Call inductor on FGC = QM0 u MM1") << flush;
    LOG(logDEBUG,*_log) << (format("  o |QM0|, |MM1|, |MM2|        %1$d %2$d %3$d") 
            % _ptop->QM0().size() % _ptop->MM1().size() % _ptop->MM2().size()).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Polarization cut-off:      ")).str() << _polar_cutoff << " nm " << flush;
    LOG(logDEBUG,*_log) << (format("  o With induction:            %1$s") % ((_polar_do_induce) ? "yes" : "no")) << flush;
    LOG(logDEBUG,*_log) << (format("  o Thole sharpness parameter: ")).str() << _polar_aDamp << flush;
    LOG(logDEBUG,*_log) << (format("  o SOR mixing factor:         ")).str() << _polar_wSOR_N << " (N) " << _polar_wSOR_C << " (C) "  << flush;
    LOG(logDEBUG,*_log) << (format("  o Iterations (max):          512")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Tolerance (rms, e*nm):     0.001")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Induce within QM0:         yes")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Subthreads:                single")).str() << flush;
    
    // return; // OVERRIDE
    
    // Forge XJob object to comply with XInductor interface
    bool polar_has_permanent_fields = true;
    XJob polar_xjob = XJob(_ptop, polar_has_permanent_fields);
    
    // INITIALIZE XINDUCTOR
    bool    polar_induce_intra_pair = true;
    int     polar_subthreads = 1;
    double  polar_epstol = 0.001;
    int     polar_maxIter = 512;
    bool    polar_maverick = _log->isMaverick(); // TODO Extract from _log
    
    XInductor polar_xind = XInductor(_polar_do_induce, 
                                     polar_induce_intra_pair, 
                                     polar_subthreads,
                                     _polar_wSOR_N,
                                     _polar_wSOR_C,
                                     polar_epstol,
                                     polar_maxIter,
                                     _polar_aDamp,
                                     polar_maverick,
                                     _top);
    polar_xind.setLog(_log);
    polar_xind.Evaluate(&polar_xjob);
    
    // SAVE RESULTS
    _polar_ETT = polar_xjob.getETOT();
    _polar_EPP = polar_xjob.getEPP();
    _polar_EPU = polar_xjob.getEPU();
    _polar_EUU = polar_xjob.getEUU();
    return;
}
        
        
void Ewald3DnD::EvaluateEnergy() {
    
    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    double EPP_fgC_mgN = ConvergeRealSpaceSum();    
    
    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    double EKK_fgC_bgP = ConvergeReciprocalSpaceSum();       
    
    // K=0 TERM (FOR 3D2D)
    double EK0_fgC_bgP = CalculateK0Correction();
    
    // SHAPE-CORRECTION (FOR 3D3D)
    double EJ_fgC_bgP = CalculateShapeCorrection();    
    
    // REAL-SPACE HIGHER-RANK CORRECTION (3D2D && 3D3D)
    double EDQ_fgC_mgN = CalculateHigherRankCorrection();
    
    // FOREGROUND CORRECTION (3D2D && 3D3D)
    double EPP_fgC_fgN = CalculateForegroundCorrection();    
    
    _ER  = EPP_fgC_mgN * _actor.int2eV;
    _EK  = EKK_fgC_bgP * _actor.int2eV;
    _E0  = EK0_fgC_bgP * _actor.int2eV;
    _EJ  = EJ_fgC_bgP  * _actor.int2eV;
    _EDQ = EDQ_fgC_mgN * _actor.int2eV;
    _EC  = EPP_fgC_fgN * _actor.int2eV;    
    _ET  = _ER + _EK + _E0 + _EJ + _EDQ - _EC;
    
    return;
}


double Ewald3DnD::ConvergeRealSpaceSum() {
    
    LOG(logDEBUG,*_log) << flush;

    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // REAL-SPACE CONVERGENCE   
    double dR = 0.1; // radius increment [nm]
    _converged_R = false;
    double prev_ER = 0.0;
    double this_ER = 0.0;
    for (int i = 0; i < 1000; ++i) {
        double Rc = _R_co + (i-1)*dR;
        // Set-up midground
        this->SetupMidground(Rc);
        // Calculate interaction energy
        this_ER = 0.;
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        this_ER += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _alpha);
                    }
                }
            }
        }
        double dER_rms = sqrt((this_ER-prev_ER)*(this_ER-prev_ER))*_actor.int2eV;
        LOG(logDEBUG,*_log)
            << (format("Rc = %1$+1.7f   |MGN| = %3$5d nm   ER = %2$+1.7f eV   dER(rms) = %4$+1.7f") 
            % Rc % (this_ER*_actor.int2eV) % _mg_N.size() % dER_rms).str() << flush;
        if (i > 0 && dER_rms < _crit_dE) {
            _converged_R = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of Rc = %1$+1.3f nm") 
                % Rc ) << flush;
            break;
        }
        prev_ER = this_ER;
    }
    return this_ER;
}


double Ewald3DnD::CalculateForegroundCorrection() {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EPP_fgC_fgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _actor.E_QQ_ERF(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    return EPP_fgC_fgN;
}


double Ewald3DnD::CalculateHigherRankCorrection() {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EDQ_fgC_mgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    EDQ_fgC_mgN += _actor.E_Q0_DQ(*(*pit1), *(*pit2));
                }
            }
        }
    }
    return EDQ_fgC_mgN;
}


string Ewald3DnD::GenerateErrorString() {
    string rstr;
    rstr += (format("Converged R-sum = %1$s, converged K-sum = %2$s")
        % ((_converged_R) ? "true" : "false")
        % ((_converged_K) ? "true" : "false")).str();
    return rstr;
}


Property Ewald3DnD::GenerateOutputString() {
    
    Property prop;
    Property &out = prop.add("output","");
    Property *next = NULL;    
    
    next = &out.add("summary", "");
    next->add("type", _jobType);
    next->add("xyz", (format("%1$+1.7f %2$+1.7f %3$+1.7f") 
        % _center.getX() % _center.getY() % _center.getZ()).str())
        .setAttribute("unit","nm");
    next->add("total", (format("%1$+1.7f") 
        % _ET).str())
        .setAttribute("unit","eV");
    
    next = &out.add("splitting", "");
    next->add("R-term", (format("%1$+1.7f") % _ER).str());
    next->add("K-term", (format("%1$+1.7f") % _EK).str());
    next->add("O-term", (format("%1$+1.7f") % _E0).str());
    next->add("J-term", (format("%1$+1.7f") % _EJ).str());
    next->add("C-term", (format("%1$+1.7f") % _EC).str());
    next->add("Q-term", (format("%1$+1.7f") % _EDQ).str());
    
    next = &out.add("shells", "");
    next->add("FGC", (format("%1$d") % _fg_C.size()).str());
    next->add("FGN", (format("%1$d") % _fg_N.size()).str());
    next->add("MGN", (format("%1$d") % _mg_N.size()).str());
    next->add("BGN", (format("%1$d") % _bg_N.size()).str());
    next->add("BGP", (format("%1$d") % _bg_P.size()).str());
    
    return prop;
}
    
    
    
    
    
    
}}