#include <votca/ctp/ewald2d.h>
#include <boost/format.hpp>


using boost::format;


namespace votca { namespace ctp {
    
  
Ewald2D::Ewald2D(Topology *top, PolarTop *ptop, XJob *xjob, Logger *log) 
    : _top(top), _ptop(ptop), _log(log), _K_co(15), _R_co(1.2) {
    
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
    
    // CALCULATE COG POSITIONS       
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
    }
    
    // SET-UP MIDGROUND
    // TODO Extend this two several molecules into the foreground
    assert(_fg_C.size() == 1);
    assert(_fg_N.size() == _fg_C.size());
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        PolarSeg *pseg = *sit;
        pseg->CalcPos();
        double dR = abs(top->PbShortestConnect(pseg->getPos(),_fg_C[0]->getPos()));
        if (dR <= _R_co) _mg_N.push_back(pseg);        
    }
    
    // CALCULATE NET CHARGES
    double Q_fg_C = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;    
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        Q_fg_C += (*sit)->CalcTotQ();
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        Q_fg_N += (*sit)->CalcTotQ();
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        Q_bg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        Q_bg_P += (*sit)->CalcTotQ();
    }
    LOG(logINFO,*_log)
        << (format("Net ground charge and size")).str()
        << flush << (format("    Q(FGC) = %1$+1.3fe |FGC| = %2$+4d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("    Q(FGN) = %1$+1.3fe |FGN| = %2$+4d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("    Q(MGN) = %1$+1.3fe |MGN| = %2$+4d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("    Q(BGN) = %1$+1.3fe |BGN| = %2$+4d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("    Q(BGP) = %1$+1.3fe |BGP| = %2$+4d") % Q_bg_P % _bg_P.size()).str()
        << flush;
    
    // COORDINATE OUTPUT FOR CHECK
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


void Ewald2D::Evaluate() {
    
    _actor.ResetEnergy();
    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    double EPP_fgC_mgN = 0.0;
    double EPP_fgC_bgN = 0.0;
    
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_mgN += _actor.E_QQ_ERFC(*(*pit1), *(*pit2));
                }                
            }           
        }
    }
    
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _bg_N.begin(); sit2 < _bg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_bgN += _actor.E_QQ_ERFC(*(*pit1), *(*pit2));
                }                
            }
        }
    }
    
    LOG(logINFO,*_log)
        << (format("Interaction FGC<-***")).str()
        << flush << (format("    EPP(FGC<-MGN) = %1$+1.7feV") % (EPP_fgC_mgN*int2eV)).str()
        << flush << (format("    EPP(FGC<-BGN) = %1$+1.7feV") % (EPP_fgC_bgN*int2eV)).str()
        << flush;
    
    return;
}
    
    
    
    
    
    
}}