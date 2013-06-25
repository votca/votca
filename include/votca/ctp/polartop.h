#ifndef __POLARTOP__H
#define	__POLARTOP__H

#include <votca/ctp/polarseg.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

    
class PolarTop
{

public:
    
   PolarTop(Topology *top) : _top(top) { cout << endl << "CONSTRUCT POLARTOP - DONE" << endl;};
  ~PolarTop();
    
   vector<PolarSeg*> &QM0() { return _qm0; }
   vector<PolarSeg*> &MM1() { return _mm1; }
   vector<PolarSeg*> &MM2() { return _mm2; }
   
   void   setQM0(vector<PolarSeg*> &qm0) { _qm0 = qm0; }
   void   setMM1(vector<PolarSeg*> &mm1) { _mm1 = mm1; }
   void   setMM2(vector<PolarSeg*> &mm2) { _mm2 = mm2; }
   
   void   CenterAround(const vec &center);
   string ShellInfoStr();
   void   PrintInfo(ostream &out);
   void   PrintPDB(string outfile);
   void   PrintInduState(string out, string format, bool split, double space);
   void   PrintInduState(FILE *out,  string format, bool split, double space);
   
    
    
private:
    
   Topology *_top; // => Box info for ::CenterAround(...)
    
   vector<PolarSeg*> _qm0; // QM SHELL
   vector<PolarSeg*> _mm1; // POLARIZABLE SHELL
   vector<PolarSeg*> _mm2; // STATIC SHELL

};  


}}

#endif





