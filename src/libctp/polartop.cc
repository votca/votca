#include <votca/ctp/polartop.h>
#include <fstream>


namespace votca { namespace ctp {

    
PolarTop::PolarTop(Topology *top) : _top(top) {
    _clean_qm0 = _clean_mm1 = _clean_mm2 = true;
    _clean_bgN = _clean_fgN = _clean_fgC = true;
};
    
    
PolarTop::~PolarTop() {
    vector<PolarSeg*> ::iterator psit;
    if (_clean_qm0) {
        for (psit = _qm0.begin(); psit < _qm0.end(); ++psit) {          
          delete *psit;          
        }
    }
    if (_clean_mm1) {
        for (psit = _mm1.begin(); psit < _mm1.end(); ++psit) {
          delete *psit;
        }
    }
    if (_clean_mm2) {
        for (psit = _mm2.begin(); psit < _mm2.end(); ++psit) {
          delete *psit;
        }
    }
    if (_clean_bgN) {
        for (psit = _bgN.begin(); psit < _bgN.end(); ++psit) {          
          delete *psit;          
        }
    }
    if (_clean_fgN) {
        for (psit = _fgN.begin(); psit < _fgN.end(); ++psit) {
          delete *psit;
        }
    }
    if (_clean_fgC) {
        for (psit = _fgC.begin(); psit < _fgC.end(); ++psit) {
          delete *psit;
        }
    }

    _qm0.clear(); _mm1.clear(); _mm2.clear();
    _bgN.clear(); _fgN.clear(); _fgC.clear();
}


void PolarTop::PrintInfo(ostream &out) {
   out  << endl << "Shells "
        << "|QM0| = " << _qm0.size() << ", "
        << "|MM1| = " << _mm1.size() << ", "
        << "|MM2| = " << _mm2.size() << ". "
        << "|FGC| = " << _fgC.size() << ", "
        << "|FGN| = " << _fgN.size() << ", "
        << "|BGN| = " << _bgN.size() << ". "
        << flush;
}


string PolarTop::ShellInfoStr() {    
    ostringstream stream;
    stream << "Shells "
           << "|QM0| = " << _qm0.size() << ", "
           << "|MM1| = " << _mm1.size() << ", "
           << "|MM2| = " << _mm2.size() << ". "
           << "|FGC| = " << _fgC.size() << ", "
           << "|FGN| = " << _fgN.size() << ", "
           << "|BGN| = " << _bgN.size() << ". ";
    return stream.str();
}


void PolarTop::Translate(const vec &shift) {
    vector<PolarSeg*> ::iterator sit;    
    for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        pseg->Translate(shift);        
    }    
    for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;     
        pseg->Translate(shift);        
    }    
    for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {        
        PolarSeg* pseg = *sit;           
        pseg->Translate(shift);
    }
    
    for (sit = _bgN.begin(); sit < _bgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        pseg->Translate(shift);        
    }    
    for (sit = _fgN.begin(); sit < _fgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;              
        pseg->Translate(shift);        
    }    
    for (sit = _fgC.begin(); sit < _fgC.end(); ++sit) {        
        PolarSeg* pseg = *sit;              
        pseg->Translate(shift);
    }
}


void PolarTop::CenterAround(const vec &center) {
    
    vector<PolarSeg*> ::iterator sit;
    
    for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);
        pseg->Translate(shift);        
    }    
    for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);        
        pseg->Translate(shift);        
    }    
    for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);        
        pseg->Translate(shift);
    }
    
    for (sit = _bgN.begin(); sit < _bgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);
        pseg->Translate(shift);        
    }    
    for (sit = _fgN.begin(); sit < _fgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);        
        pseg->Translate(shift);        
    }    
    for (sit = _fgC.begin(); sit < _fgC.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        vec shift = _top->PbShortestConnect(center, pseg->getPos())
                         -(pseg->getPos() - center);        
        pseg->Translate(shift);
    }
    
    _center = center;
}


void PolarTop::PrintPDB(string outfile) {
    
    FILE *out;
    out = fopen(outfile.c_str(),"w");
    vector<PolarSeg*> ::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "QM0");
        }
    }
    for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM1");
        }
    }
    for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM2");
        }
    }
    for (sit = _fgC.begin(); sit < _fgC.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGC");
        }
    }
    for (sit = _fgN.begin(); sit < _fgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGN");
        }
    }
    for (sit = _bgN.begin(); sit < _bgN.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGN");
        }
    }   
    fclose(out);    
}


void PolarTop::PrintInduState(string outfile, string format,
                              bool split_dpl, double dpl_space) {
    
    FILE *out;
    out = fopen(outfile.c_str(), "w");
    this->PrintInduState(out, format, split_dpl, dpl_space);
    fclose(out);    
}


void PolarTop::PrintInduState(FILE *out, string format,
                              bool split_dpl, double dpl_space) {

    vector< PolarSeg* > ::iterator sit;
    vector< APolarSite* > ::iterator pit;    

    // Count polar sites for header line in xyz
    int pcount = 0;
    if (format == "xyz") {
        for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {
            pcount += (*sit)->size();
        }
        for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {
            int mult = (split_dpl) ? 3 : 1;
            pcount += mult*(*sit)->size();
        }
        for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {
            pcount += (*sit)->size();
        }

        fprintf(out, "%1d\nXYZ WITH DIPOLES SPLIT FOR SHELL MM1\n", pcount);
    }

    
    // COORDINATES OF QM0 SHELL
    for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {
        vec pb_shift = vec(0,0,0);
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->WriteXyzLine(out, pb_shift, format);
        }
    }
    
    if (format == "gaussian") fprintf(out, "\n");
    
    // INDUCTION STATE OF MM1
    for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {
        vec pb_shift = vec(0,0,0);
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->WriteChkLine(out, pb_shift, split_dpl, format, dpl_space);
        }
    }

    // POINT CHARGES OF MM2
    for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {
        vec pb_shift = vec(0,0,0);
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->WriteChkLine(out, pb_shift, false, format, dpl_space);
        }
    }
    
}

}}
