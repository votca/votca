#include <votca/ctp/polartop.h>
#include <fstream>


namespace votca { namespace ctp {

PolarTop::~PolarTop() {
       
      vector<PolarSeg*> ::iterator psit;
      for (psit = _qm0.begin(); psit < _qm0.end(); ++psit) {          
          delete *psit;          
      }
      for (psit = _mm1.begin(); psit < _mm1.end(); ++psit) {
          delete *psit;
      }
      for (psit = _mm2.begin(); psit < _mm2.end(); ++psit) {
          delete *psit;
      }
      
      _qm0.clear(); _mm1.clear(); _mm2.clear();       
      
      cout << endl << "Destructed ptop" << endl;
      
   }    

void PolarTop::PrintInfo() {

   cout << endl
        << "|QM0| = " << _qm0.size() << ", "
        << "|MM1| = " << _mm1.size() << ", "
        << "|MM2| = " << _mm2.size() << ". "
        << flush;
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
}


void PolarTop::PrintPDB(const string &outfile) {
    
    FILE *out;
    out = fopen(outfile.c_str(),"w");
    vector<PolarSeg*> ::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _qm0.begin(); sit < _qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->PSites().begin(); pit < pseg->PSites().end(); ++pit) {
            (*pit)->WritePdbLine(out);
        }
    }
    for (sit = _mm1.begin(); sit < _mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->PSites().begin(); pit < pseg->PSites().end(); ++pit) {
            (*pit)->WritePdbLine(out);
        }
    }
    for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->PSites().begin(); pit < pseg->PSites().end(); ++pit) {
            (*pit)->WritePdbLine(out);
        }
    }
    fclose(out);
    
}
    
}}
