#include <votca/ctp/ewald2d.h>
#include <boost/format.hpp>
#include <algorithm>


using boost::format;


namespace votca { namespace ctp {
    
  
Ewald2D::Ewald2D(Topology *top, PolarTop *ptop, XJob *xjob, double R_co, Logger *log) 
    : _top(top), _ptop(ptop), _log(log) {
    
    // EWALD INTERACTION PARAMETERS
    _R_co = R_co;
    _K_co = 3.5*3.5/_R_co;
    _a = 3.5/_R_co;
    
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
    
    
    // SET-UP MIDGROUND
    // TODO Extend this to several molecules in the foreground
    assert(_fg_C.size() == 1);
    assert(_fg_N.size() == _fg_C.size());
    Segment *central = _top->getSegment(_fg_C[0]->getId());
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        PolarSeg *pseg = *sit;
        Segment *other = _top->getSegment(pseg->getId());
        double dR = abs(top->PbShortestConnect(other->getPos(),central->getPos()));
        if (dR <= _R_co) _mg_N.push_back(pseg);
    }
    
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
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Net ground charge and size")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+4d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+4d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| = %2$+4d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+4d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+4d") % Q_bg_P % _bg_P.size()).str()
        << flush;
    
    // CHARGE NEUTRAL & DEPOLARIZE, COORDINATE OUTPUT FOR CHECK
    FILE *out;
    string outfile = xjob->getTag()+"_FGC_FGN_MGN_BGN_BGP.pdb";
    out = fopen(outfile.c_str(),"w");
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGC");
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGN");
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MGN");
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGN");
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGP");
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }
    }
    fclose(out);    
    
    // INIT. INTERACTOR
    _actor = XInteractor(top, 0.39);
    
    return;
}



void Ewald2D::CheckParameters() {
    
    LOG(logDEBUG,*_log) << flush;

    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // REAL-SPACE CONVERGENCE
    vector< PolarSeg* > mg;    
    double dR = 0.1; // [nm]
    for (int i = 0; i < 21; ++i) {
        double Rc = _R_co + (i-10)*dR;
        // Set-up midground
        mg.clear();
        Segment *centralSeg = _top->getSegment(_fg_C[0]->getId());
        for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
            PolarSeg *pseg = *sit;
            Segment *other = _top->getSegment(pseg->getId());
            pseg->CalcPos();
            double dR = abs(_top->PbShortestConnect(other->getPos(),
                centralSeg->getPos()));
            if (dR <= Rc) mg.push_back(pseg);
        }
        // Reset interactor
        _actor.ResetEnergy();
        double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
        // Calculate interaction energy
        double epp = 0.;
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = mg.begin(); sit2 < mg.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        epp += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _a);
                    }
                }
            }
        }
        LOG(logINFO,*_log)
            << (format("ERR(Rc=%1$+1.7fnm, |.|=%3$04d) = %2$+1.7f eV") 
            % Rc % (epp*int2eV) % mg.size()).str() << flush;
    }  
    
}


void Ewald2D::Evaluate() {
    
    _actor.ResetEnergy();
    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;    
    
    // SET-UP RECIPROCAL LATTICE
    vec a = _top->getBox().getCol(0);
    vec b = _top->getBox().getCol(1);
    vec c = _top->getBox().getCol(2);    
    double LxLyLz = a*(b^c);
    double LxLy = abs(a ^ b);
    
    vec A = 2*M_PI/LxLyLz * b^c;
    vec B = 2*M_PI/LxLyLz * c^a;
    vec C = 2*M_PI/LxLyLz * a^b;
    double KA = abs(A);
    double KB = abs(B);
    double KC = abs(C);
    
    // Careful, this is not quite ideal for non-orthogonal A & B
    int kx_max = ceil(_K_co/KA);
    int ky_max = ceil(_K_co/KB);
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "System & Ewald parameters" << flush;
    LOG(logDEBUG,*_log) << "  o Real-space unit cell:    " << a << b << c << flush;
    LOG(logDEBUG,*_log) << "  o Real-space c/o:          " << _R_co << " nm" << flush;
    LOG(logDEBUG,*_log) << "  o 1st Brillouin zone:      " << A << B << C << flush;
    LOG(logDEBUG,*_log) << "  o Reciprocal-space c/o:    " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o R-K switching param.     " << _a << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o Unit-cell volume:        " << LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*_log) << "  o LxLy for 2D Ewald:       " << LxLy << " nm**2" << flush;  
    LOG(logDEBUG,*_log) << "  o kx(max), ky(max):        " << kx_max << ", " << ky_max << flush;
    
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
    int kx = 0;
    int kz = 0; // This is 2D Ewald
    int ky = 0;
    for (ky = 1; ky < ky_max+1; ++ky) {        
        vec K = kx*A + ky*B;
        if (abs(K)>_K_co) { continue; }
        Ks.push_back(K);
    }
    for (kx = 1; kx < kx_max+1; ++kx) {
        for (ky = -ky_max; ky < ky_max+1; ++ky) {
            vec K = kx*A + ky*B;            
            if (abs(K)>_K_co) { continue; }            
            Ks.push_back(K);
        }
    }
    std::sort(Ks.begin(), Ks.end(), vec_smaller_than());
    
    // REAL-SPACE CONTRIBUTION
    double EPP_fgC_mgN = 0.0;
    double EPP_fgC_bgN = 0.0;    
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_mgN += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _a);
                }                
            }           
        }
    }    
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _bg_N.begin(); sit2 < _bg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_bgN += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _a);
                }                
            }
        }
    }
    
    // K=0 TERM
    double EK0_fgC_bgP = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EK0_fgC_bgP += _actor.E_QQ_K0(*(*pit1), *(*pit2), _a);
                }
            }
        }
    }
    EK0_fgC_bgP *= 2*M_PI/LxLy;
    
    // K=K TERM
    LOG(logDEBUG,*_log) << flush;
    double EKK_fgC_bgP = 0.0;    
    int N_EKK_memory = int(0.5*(kx_max+ky_max)+0.5);
    int N_K_proc = 0;
    vector< double > dEKKs;
    double dEKK_rms_crit = 1e-3;
    
    for (kit = Ks.begin(); kit < Ks.end(); ++kit, ++N_K_proc) {
        vec k = *kit;
        double K = abs(k);
        double dEKK = 0.0;
        LOG(logDEBUG,*_log)  
            << (format("Cell = %1$+1.1f %2$+1.1f %3$+1.1f   |K| = %4$+1.3f 1/nm") 
            % (k.getX()/KA) % (k.getY()/KB) % (k.getZ()/KC) % K);        
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {                                                
                        dEKK += _actor.E_QQ_KK(*(*pit1), *(*pit2), _a, k);
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
            << (format("   dE = %1$+1.7f   EKK = %2$+1.7f   RMS(%4$d) = %3$+1.7f") % dEKK % EKK_fgC_bgP % dEKK_rms % N_EKK_memory) << flush;
        
        if (dEKK_rms <= 1e-4) {
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged as of |K| >= %1$+1.3f 1/nm") % K ) << flush;
            break;
        }
        
        
    }
    EKK_fgC_bgP *= 2*M_PI/LxLy;
    
    // FOREGROUND CORRECTION
    double EPP_fgC_fgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _actor.E_QQ_ERF(*(*pit1), *(*pit2), _a);
                }
            }
        }
    }
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  o EPP(FGC->BGN) = %1$+1.7f eV") % (EPP_fgC_bgN*int2eV)).str()
        << flush << (format("  + EPP(FGC->MGN) = %1$+1.7f eV") % (EPP_fgC_mgN*int2eV)).str()
        << flush << (format("  + EK0(FGC->BGP) = %1$+1.7f eV") % (EK0_fgC_bgP*int2eV)).str()
        << flush << (format("  + EKK(FGC->BGP) = %1$+1.7f eV") % (EKK_fgC_bgP*int2eV)).str()
        << flush << (format("  - EPP(FGC->FGN) = %1$+1.7f eV") % (EPP_fgC_fgN*int2eV)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)        = %1$+1.7f eV") % ((EPP_fgC_mgN+EK0_fgC_bgP+EKK_fgC_bgP-EPP_fgC_fgN)*int2eV)).str()
        << flush;
    LOG(logDEBUG,*_log) << flush;
    
    return;
}
    
    
    
    
    
    
}}