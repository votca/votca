#ifndef __POLARTOP__H
#define	__POLARTOP__H

#include <votca/xtp/polarseg.h>
#include <votca/xtp/topology.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace votca { namespace xtp {

    
class PolarTop
{

public:
   
   PolarTop();
   PolarTop(Topology *top);
  ~PolarTop();
   
   // QM-MM-MM
   vector<PolarSeg*> &QM0() { return _qm0; }
   vector<PolarSeg*> &MM1() { return _mm1; }
   vector<PolarSeg*> &MM2() { return _mm2; }   
   
   void   setQM0(vector<PolarSeg*> &qm0, bool clean = true) { _qm0 = qm0; _clean_qm0 = clean; }
   void   setMM1(vector<PolarSeg*> &mm1, bool clean = true) { _mm1 = mm1; _clean_mm1 = clean; }
   void   setMM2(vector<PolarSeg*> &mm2, bool clean = true) { _mm2 = mm2; _clean_mm2 = clean; }
   
   void   setSegsQM0(vector<Segment*> &qm0) { _segs_qm0 = qm0; }
   void   setSegsMM1(vector<Segment*> &mm1) { _segs_mm1 = mm1; }
   void   setSegsMM2(vector<Segment*> &mm2) { _segs_mm2 = mm2; }
   
   // EWALD FGC-FGN-BGN
   vector<PolarSeg*> &BGN() { return _bgN; }
   vector<PolarSeg*> &FGN() { return _fgN; }
   vector<PolarSeg*> &FGC() { return _fgC; }
   
   void   setBGN(vector<PolarSeg*> &bgN, bool clean = true) { _bgN = bgN; _clean_bgN = clean; }
   void   setFGN(vector<PolarSeg*> &fgN, bool clean = true) { _fgN = fgN; _clean_fgN = clean; }
   void   setFGC(vector<PolarSeg*> &fgC, bool clean = true) { _fgC = fgC; _clean_fgC = clean; }
   
   void   setSegsBGN(vector<Segment*> &bgN) { _segs_bgN = bgN; }
   void   setSegsFGN(vector<Segment*> &fgN) { _segs_fgN = fgN; }
   void   setSegsFGC(vector<Segment*> &fgC) { _segs_fgC = fgC; }
   
   // TRANSFORMATIONS & POINT OF REFERENCE
   void   Translate(const vec &shift);
   void   CenterAround(const vec &center);
   vec    getCenter() { return _center; }
   
   string ShellInfoStr();
   void   PrintInfo(ostream &out);
   void   PrintPDB(string outfile);
   void   PrintInduState(string out, string format, bool split, double space);
   void   PrintInduState(FILE *out,  string format, bool split, double space);
   
   // POLARIZATION STATE
   void   setPolarizationIter(int iter, bool converged) { 
       _polarization_iter = iter; _polarization_converged = converged; }
   int    getPolarizationIter() { return _polarization_iter; }
   
   template<class Archive>
   void serialize(Archive &arch, const unsigned int version) {
       arch & _center;       
       arch & _qm0;
       arch & _mm1;
       arch & _mm2;        
       arch & _bgN;
       arch & _fgN;
       arch & _fgC;
       if (version > 0) {
           arch & _polarization_iter;
           arch & _polarization_converged;
       }
       return;
   }
   
    void RemoveAllOwnership();
    void LoadFromDrive(string archfile);
    void SaveToDrive(string archfile);
    
private:
    
   Topology *_top; // => Box info for ::CenterAround(...)
   vec _center;    // => Set in ::CenterAround(...)
   
   // QM-MM-MM
   vector<PolarSeg*> _qm0; // QM SHELL
   vector<PolarSeg*> _mm1; // POLARIZABLE SHELL
   vector<PolarSeg*> _mm2; // STATIC SHELL
   
   vector<Segment*> _segs_qm0;
   vector<Segment*> _segs_mm1;
   vector<Segment*> _segs_mm2;
   
   bool _clean_qm0;
   bool _clean_mm1;
   bool _clean_mm2;
   
   // EWALD FGC-FGN-BGN
   vector<PolarSeg*> _bgN; // Neutral background
   vector<PolarSeg*> _fgN; // Neutral foreground
   vector<PolarSeg*> _fgC; // Charged foreground
   
   vector<Segment*> _segs_bgN;
   vector<Segment*> _segs_fgN;
   vector<Segment*> _segs_fgC;
   
   bool _clean_bgN;
   bool _clean_fgN;
   bool _clean_fgC;
   
   int _polarization_iter;
   bool _polarization_converged;

};

}}

#endif





