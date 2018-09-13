/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9


#ifndef VOTCA_XTP_POLARTOP_H
#define	VOTVA_XTP_POLARTOP_H

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
   std::vector<PolarSeg*> &QM0() { return _qm0; }
   std::vector<PolarSeg*> &MM1() { return _mm1; }
   std::vector<PolarSeg*> &MM2() { return _mm2; }   
   
   void   setQM0(std::vector<PolarSeg*> &qm0, bool clean = true) { _qm0 = qm0; _clean_qm0 = clean; }
   void   setMM1(std::vector<PolarSeg*> &mm1, bool clean = true) { _mm1 = mm1; _clean_mm1 = clean; }
   void   setMM2(std::vector<PolarSeg*> &mm2, bool clean = true) { _mm2 = mm2; _clean_mm2 = clean; }
   
   void   setSegsQM0(std::vector<Segment*> &qm0) { _segs_qm0 = qm0; }
   void   setSegsMM1(std::vector<Segment*> &mm1) { _segs_mm1 = mm1; }
   void   setSegsMM2(std::vector<Segment*> &mm2) { _segs_mm2 = mm2; }
   
   // EWALD FGC-FGN-BGN
   std::vector<PolarSeg*> &BGN() { return _bgN; }
   std::vector<PolarSeg*> &FGN() { return _fgN; }
   std::vector<PolarSeg*> &FGC() { return _fgC; }
   
   void   setBGN(std::vector<PolarSeg*> &bgN, bool clean = true) { _bgN = bgN; _clean_bgN = clean; }
   void   setFGN(std::vector<PolarSeg*> &fgN, bool clean = true) { _fgN = fgN; _clean_fgN = clean; }
   void   setFGC(std::vector<PolarSeg*> &fgC, bool clean = true) { _fgC = fgC; _clean_fgC = clean; }
   
   void   setSegsBGN(std::vector<Segment*> &bgN) { _segs_bgN = bgN; }
   void   setSegsFGN(std::vector<Segment*> &fgN) { _segs_fgN = fgN; }
   void   setSegsFGC(std::vector<Segment*> &fgC) { _segs_fgC = fgC; }
   
   // TRANSFORMATIONS & POINT OF REFERENCE
   void   Translate(const vec &shift);
   void   CenterAround(const vec &center);
   vec    getCenter() { return _center; }
   
   std::string ShellInfoStr();
   void   PrintInfo(std::ostream &out);
   void   PrintPDB(std::string outfile);
   void   PrintInduState(std::string out, std::string format, bool split, double space);
   void   PrintInduState(FILE *out,  std::string format, bool split, double space);
   
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
    void LoadFromDrive(std::string archfile);
    void SaveToDrive(std::string archfile);
    
private:
    
   Topology *_top; // => Box info for ::CenterAround(...)
   vec _center;    // => Set in ::CenterAround(...)
   
   // QM-MM-MM
   std::vector<PolarSeg*> _qm0; // QM SHELL
   std::vector<PolarSeg*> _mm1; // POLARIZABLE SHELL
   std::vector<PolarSeg*> _mm2; // STATIC SHELL
   
   std::vector<Segment*> _segs_qm0;
   std::vector<Segment*> _segs_mm1;
   std::vector<Segment*> _segs_mm2;
   
   bool _clean_qm0;
   bool _clean_mm1;
   bool _clean_mm2;
   
   // EWALD FGC-FGN-BGN
   std::vector<PolarSeg*> _bgN; // Neutral background
   std::vector<PolarSeg*> _fgN; // Neutral foreground
   std::vector<PolarSeg*> _fgC; // Charged foreground
   
   std::vector<Segment*> _segs_bgN;
   std::vector<Segment*> _segs_fgN;
   std::vector<Segment*> _segs_fgC;
   
   bool _clean_bgN;
   bool _clean_fgN;
   bool _clean_fgC;
   
   int _polarization_iter;
   bool _polarization_converged;

};

}}

#endif // VOTCA_XTP_POLARTOP_H
