#include <votca/ctp/ewald2d.h>
#include <boost/format.hpp>
#include <algorithm>


using boost::format;


namespace votca { namespace ctp {


Ewald2D::~Ewald2D() {
    vector< PolarSeg* >::iterator sit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
}
    
    
Ewald2D::Ewald2D(Topology *top, PolarTop *ptop, double R_co, Logger *log) 
    : _top(top), _ptop(ptop), _log(log) {
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _R_co = R_co;  
    _K_co = 100/_R_co;
    _alpha = 3.5/_R_co;
    
    
    // SET-UP REAL & RECIPROCAL SPACES
    _a = _top->getBox().getCol(0);
    _b = _top->getBox().getCol(1);
    _c = _top->getBox().getCol(2);
    _LxLyLz = _a*(_b^_c);
    _LxLy = abs(_a ^ _b);
    
    _A = 2*M_PI/_LxLyLz * _b^_c;
    _B = 2*M_PI/_LxLyLz * _c^_a;
    _C = 2*M_PI/_LxLyLz * _a^_b;
    double KA = abs(_A);
    double KB = abs(_B);
    double KC = abs(_C);
    
    // .ff: This may not be sufficient for non-orthogonal a,b
    _na_max = ceil(_R_co/abs(_a)-0.5)+1;
    _nb_max = ceil(_R_co/abs(_b)-0.5)+1;
    
    // .ff: This may not be sufficient for non-orthogonal A,B
    _NA_max = ceil(_K_co/KA);
    _NB_max = ceil(_K_co/KB); 
    
    
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
    
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;  
    
    // SET-UP MIDGROUND (INCLUDING PERIODIC IMAGES IF REQUIRED)
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Generate periodic images. ";
    this->SetupMidground(_R_co);
    
    // CALCULATE COG POSITIONS, NET CHARGE
    double Q_fg_C = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_C += (*sit)->CalcTotQ();
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_N += (*sit)->CalcTotQ();
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }   
    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| = %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush;
    
    
    // CHARGE APPROPRIATELY & DEPOLARIZE
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    
    // CALCULATE NET DIPOLE OF BGP
    vec netdpl = vec(0,0,0);
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl += (*pit)->getPos() * (*pit)->getQ00();
        }
    }
    LOG(logINFO,*_log)
        << (format("Net dipole moment of background density")).str()
        << flush << (format("  o D(BGP) [e*nm] = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl.getX() % netdpl.getY() % netdpl.getZ()).str() << flush;
    
    // INIT. INTERACTOR
    _actor = XInteractor(top, 0.39);    
    return;
}


void Ewald2D::SetupMidground(double R_co) {    
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
    for (int na = -_na_max; na < _na_max+1; ++na)
        nas.push_back(na);
    for (int nb = -_nb_max; nb < _nb_max+1; ++nb)
        nbs.push_back(nb);
    vector< int >::iterator ait;
    vector< int >::iterator bit;
    
    // SAMPLE MIDGROUND FROM BGP EXCLUDING CENTRAL SEG.
    assert(_fg_N.size() == _fg_C.size());
    Segment *central = _top->getSegment(_fg_C[0]->getId());
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        Segment *other = _top->getSegment(pseg->getId());        
        // Periodic images
        for (ait = nas.begin(); ait < nas.end(); ++ait) {
            int na = *ait;
            for (bit = nbs.begin(); bit < nbs.end(); ++bit) {                
                int nb = *bit;                
                if (na == 0 && nb == 0 && central->getId() == other->getId())
                    continue;
                vec L = na*_a + nb*_b;
                vec pos_L = pseg->getPos() + L;
                double dR_L = abs(central->getPos()-pos_L);
                if (dR_L <= R_co) {
                    PolarSeg *newSeg = new PolarSeg(pseg);
                    newSeg->Translate(L);
                    _mg_N.push_back(newSeg);
                }
            }
        }
    }
    return;
}


void Ewald2D::WriteDensitiesPDB(string pdbfile) {
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


double Ewald2D::ConvergeRealSpaceSum() {
    
    LOG(logDEBUG,*_log) << flush;

    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // REAL-SPACE CONVERGENCE   
    double dR = 0.1; // [nm]
    double dE_tol = 1e-4; // [eV]
    _converged_R = false;
    double prev_ER = 0.0;
    double this_ER = 0.0;
    for (int i = 0; i < 100; ++i) {
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
        if (i > 0 && dER_rms < dE_tol) {
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


double Ewald2D::ConvergeReciprocalSpaceSum() {
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;    
    
    // CELLS OF THE RECIPROCAL LATTICE OVER WHICH TO SUM
    //    ---------------------------------
    //      | x | x | x | x |   |   |   |
    //    ---------------------------------            
    //      | x | x | x | x |   |   |   |
    //    ---------------------------------
    //      | x | x | x | O |   |   |   |
    //    ---------------------------------
    //      | x | x | x |   |   |   |   |  
    //    ---------------------------------
    //      | x | x | x |   |   |   |   |  
    //    ---------------------------------
    
    vector< vec > Ks;
    vector< vec >::iterator kit;
    int kz = 0; // This is 2D Ewald
    int kx = 0;
    int ky = 0;
    for (ky = 1; ky < _NB_max+1; ++ky) {        
        vec K = kx*_A + ky*_B;
        if (abs(K)>_K_co) { continue; }
        Ks.push_back(K);
    }
    for (kx = 1; kx < _NA_max+1; ++kx) {
        for (ky = -_NB_max; ky < _NB_max+1; ++ky) {
            vec K = kx*_A + ky*_B;            
            if (abs(K)>_K_co) { continue; }            
            Ks.push_back(K);
        }
    }
    std::sort(Ks.begin(), Ks.end(), VecSmallerThan());  
        
    // K=K TERM
    LOG(logDEBUG,*_log) << flush;
    double EKK_fgC_bgP = 0.0;    
    int N_EKK_memory = int(0.5*(_NA_max+_NB_max)+0.5);
    int N_K_proc = 0;
    vector< double > dEKKs;
    double dEKK_rms_crit = 1e-4; // [eV]
    _converged_K = false;
    
    for (kit = Ks.begin(); kit < Ks.end(); ++kit, ++N_K_proc) {
        vec k = *kit;
        double K = abs(k);
        double dEKK = 0.0;
        LOG(logDEBUG,*_log)  
            << (format("k = %1$+1.3f %2$+1.3f %3$+1.3f   |K| = %4$+1.3f 1/nm") 
            % (k.getX()) % (k.getY()) % (k.getZ()) % K);
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {                                                
                        dEKK += _actor.E_QQ_KK(*(*pit1), *(*pit2), _alpha, k);
                    }
                }
            }
        }
        
        EKK_fgC_bgP += dEKK;
        
        // CONVERGED?
        double dEKK_rms = 0.0;
        if (dEKKs.size() < N_EKK_memory) {
            dEKKs.resize(N_EKK_memory,dEKK);
        }
        else {
            dEKKs[N_K_proc % N_EKK_memory] = dEKK;
        }
        for (int i = 0; i < dEKKs.size(); ++i) {
            dEKK_rms += dEKKs[i]*dEKKs[i];
        }
        dEKK_rms /= dEKKs.size();
        dEKK_rms = sqrt(dEKK_rms);
        
        LOG(logDEBUG,*_log)  
            << (format("   dE = %1$+1.7f   EKK = %2$+1.7f   RMS(%4$d) = %3$+1.7f") 
            % (dEKK*2*M_PI/_LxLy*_actor.int2eV) % (EKK_fgC_bgP*2*M_PI/_LxLy*_actor.int2eV)
            % (dEKK_rms*2*M_PI/_LxLy*_actor.int2eV) % N_EKK_memory) << flush;
        
//        cout << endl;
//        for (int i = 0; i < dEKKs.size(); ++i) {
//            cout << "| " << dEKKs[i];
//        }
//        cout << flush;
        
        if (dEKK_rms*2*M_PI/_LxLy*_actor.int2eV <= dEKK_rms_crit && N_K_proc > 2) {
        //if (dEKK_rms*2*M_PI/LxLy*_actor.int2eV <= dEKK_rms_crit) {
            _converged_K = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of |K| = %1$+1.3f 1/nm") 
                % K ) << flush;
            break;
        }
        
        
    }
    EKK_fgC_bgP *= 2*M_PI/_LxLy;
    return EKK_fgC_bgP;    
}


void Ewald2D::Evaluate() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "System & Ewald parameters (3D x 2D)" << flush;
    LOG(logDEBUG,*_log) << "  o Real-space unit cell:      " << _a << " x " << _b << " x " << _c << flush;
    LOG(logDEBUG,*_log) << "  o Real-space c/o (guess):    " << _R_co << " nm" << flush;
    LOG(logDEBUG,*_log) << "  o na(max), nb(max):          " << _na_max << ", " << _nb_max << flush;
    LOG(logDEBUG,*_log) << "  o 1st Brillouin zone:        " << _A << " x " << _B << " x " << _C << flush;
    LOG(logDEBUG,*_log) << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o R-K switching param.       " << _alpha << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*_log) << "  o LxLy for 2D Ewald:         " << _LxLy << " nm**2" << flush;  
    LOG(logDEBUG,*_log) << "  o kx(max), ky(max):          " << _NA_max << ", " << _NB_max << flush;
        
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;    
    
    // REAL-SPACE CONTRIBUTION
    double EPP_fgC_mgN = ConvergeRealSpaceSum();    
    
    // RECIPROCAL-SPACE CONTRIBUTION
    double EKK_fgC_bgP = ConvergeReciprocalSpaceSum();       
    
    // K=0 TERM
    double EK0_fgC_bgP = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EK0_fgC_bgP += _actor.E_QQ_K0(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    EK0_fgC_bgP *= 2*M_PI/_LxLy;
    
    // FOREGROUND CORRECTION
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
    
    // REAL-SPACE HIGHER-RANK CORRECTION
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
    
    _ER = EPP_fgC_mgN*_actor.int2eV;
    _EC = EPP_fgC_fgN*_actor.int2eV;
    _EK = EKK_fgC_bgP*_actor.int2eV;
    _E0 = EK0_fgC_bgP*_actor.int2eV;
    _ET = _ER + _EK + _E0 - _EC;
    _EDQ = EDQ_fgC_mgN*_actor.int2eV;
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN) = %1$+1.7f eV") % _ER).str()
        << flush << (format("  + EK0(FGC->BGP) = %1$+1.7f eV") % _E0).str()
        << flush << (format("  + EKK(FGC->BGP) = %1$+1.7f eV") % _EK).str()
        << flush << (format("  - EPP(FGC->FGN) = %1$+1.7f eV") % _EC).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)        = %1$+1.7f eV") % _ET).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + EDQ(FGC->MGN) = %1$+1.7f eV") % _EDQ).str()
        << flush;
    LOG(logDEBUG,*_log) << flush;    
    return;
}


string Ewald2D::GenerateErrorString() {
    string rstr;
    rstr += (format("Converged R-sum = %1$s, converged K-sum = %2$s")
        % ((_converged_R) ? "true" : "false")
        % ((_converged_K) ? "true" : "false")).str();
    return rstr;
}


string Ewald2D::GenerateOutputString() {
    string rstr;
    rstr += (format("XYZ %1$+1.7f %2$+1.7f %3$+1.7f ") 
        % _center.getX() % _center.getY() % _center.getZ()).str();
    rstr += (format("ET %1$+1.7f ER %2$+1.7f EK %3$+1.7f E0 %4$+1.7f EC %5$+1.7f EDQ %6$+1.7f ") 
        % _ET % _ER % _EK % _E0 % _EC %_EDQ).str();
    rstr += (format("FGC %1$1d FGN %2$1d MGN %3$3d BGN %4$4d BGP %5$4d") 
        % _fg_C.size() % _fg_N.size() % _mg_N.size() % _bg_N.size() 
        % _bg_P.size()).str();
    return rstr;
}
    
    
    
    
    
    
}}