#include <votca/ctp/ewaldnd.h>
#include <boost/format.hpp>
#include <algorithm>


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
}
    
    
Ewald3DnD::Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
    : _top(top), _ptop(ptop), _log(log) {
    
    // EVALUATE OPTIONS
    _R_co = opt->get("options.ewald.coulombmethod.cutoff").as<double>();
    _crit_dE = opt->get("options.ewald.convergence.energy").as<double>();
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _K_co = 100/_R_co;
    _alpha = 3.5/_R_co;
    
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
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();

    _fg_C = ptop->FGC();
    _fg_N = ptop->FGN();
    _bg_N = ptop->BGN();        
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
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_C += (*sit)->CalcTotQ();
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
    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| ~ %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
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


void Ewald3DnD::SetupMidground(double R_co) {
    // SET-UP MIDGROUND
    // TODO Extend this to several molecules in the foreground
    assert(_fg_C.size() == 1);
    _center = _fg_C[0]->getPos();
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
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN)  = %1$+1.7f eV") % _ER).str()
        << flush << (format("  + EKK(FGC->BGP)  = %1$+1.7f eV") % _EK).str()       
        << flush << (format("  - EPP(FGC->FGN)  = %1$+1.7f eV") % _EC).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)         = %1$+1.7f eV") % (_ET-_EJ-_EDQ-_E0)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + EK0(FGC->BGP)  = %1$+1.7f eV") % _E0).str() 
        << flush << (format("  + EDQ(FGC->MGN)  = %1$+1.7f eV") % _EDQ).str()
        << flush << (format("  + EJ(shape-dep.) = %1$+1.7f eV") % _EJ).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + SUM(E) (0,Q,J) = %1$+1.7f eV") % (_ET)).str()
        << flush;
    LOG(logDEBUG,*_log) << flush;
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


string Ewald3DnD::GenerateOutputString() {
    string rstr;
    rstr += (format("XYZ %1$+1.7f %2$+1.7f %3$+1.7f ") 
        % _center.getX() % _center.getY() % _center.getZ()).str();
    rstr += (format("ET %1$+1.7f ER %2$+1.7f EK %3$+1.7f E0 %4$+1.7f EC %5$+1.7f EJ %6$+1.7f EDQ %7$+1.7f ") 
        % _ET % _ER % _EK % _E0 % _EC % _EJ %_EDQ).str();
    rstr += (format("FGC %1$1d FGN %2$1d MGN %3$3d BGN %4$4d BGP %5$4d") 
        % _fg_C.size() % _fg_N.size() % _mg_N.size() % _bg_N.size() 
        % _bg_P.size()).str();
    return rstr;
}
    
    
    
    
    
    
}}