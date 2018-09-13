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

#ifndef VOTCA_XTP_XJOB_H
#define	VOTCA_XTP_XJOB_H

#include <votca/xtp/topology.h>
#include <votca/xtp/job.h>

namespace votca { namespace xtp {

class PolarTop;

class XJob
{
public:

    XJob(int id, std::string tag, std::vector<Segment*> &qmSegs, 
         std::vector<std::string> &qmSegMps, Topology *top);
    
    XJob(PolarTop *ptop, bool start_from_cpt);

   ~XJob();

   int                  getId()                     { return _id; }
   std::string               getTag()                    { return _tag; }
   std::vector<Segment*>    &getSegments()               { return _qmSegs; }
   std::vector<std::string>      &getSegMps()                 { return _qmSegMps; }
   Topology            *getTop()                    { return _top; }
   PolarTop            *getPolarTop()               { return _ptop; }
   int                  getIter()                   { return _iter; }
   int                  getUserId() { return _userId; }
   void                 setUserId(int userId) { _userId = userId; }
   void                 setPolarTop(PolarTop *ptop);
   void                 setInfoLine(bool printMM = true, bool printQM = true);
   std::string               getInfoLine() { return _infoLine; }
   
   void                 CalcCenterPos();
   vec                 &Center()                    { return _center; }
   
   inline bool          isInCenter(int segId);
   inline bool          isWithinDist(const vec &pt, double dist, Topology *top);
   bool                 StartFromCPT()  { return _start_from_cpt; }
   void                 WriteInfoLine(FILE *out);
   Property             GenerateOutputProperty();

   
   
   double getEF00() { return _EF_PAIR_PAIR; }
   double getEF01() { return _EF_PAIR_SPH1; }
   double getEF02() { return _EF_PAIR_SPH2; }
   double getEF11() { return _EF_SPH1_SPH1; }
   double getEF12() { return _EF_SPH1_SPH2; }
   double getEM0()  { return _EM_PAIR; }
   double getEM1()  { return _EM_SPH1; }
   double getEM2()  { return _EM_SPH2; }
   double getETOT() { return _E_Tot; }
   
   double getEQM()  { return _E_QM; }
   double getESF()  { return _E_SF; }
   
   double getEPP()  { return _EPP; }
   double getEPU()  { return _EPU; }
   double getEUU()  { return _EUU; }
   
   void setInduIter(int iter) { 
       _iter = iter;
   }
   
   void setShellSizes(int qm0_size, int mm1_size, int mm2_size) {
       _qm0_size = qm0_size;
       _mm1_size = mm1_size;
       _mm2_size = mm2_size;
   }
   
   void     setEnergy(double E_Tot,   
                      double E_Pair_Pair, 
                      double E_Pair_Sph1, 
                      double E_Pair_Sph2,
                      double E_Sph1_Sph1,
                      double E_Sph1_Sph2,
                      double E_P,   
                      double E_U) {

       _E_Tot       = E_Tot;
       _E_Pair_Pair = E_Pair_Pair;
       _E_Pair_Sph1 = E_Pair_Sph1;
       _E_Pair_Sph2 = E_Pair_Sph2;
       _E_Sph1_Sph1 = E_Sph1_Sph1;
       _E_Sph1_Sph2 = E_Sph1_Sph2;

       _E_PERM      = E_P;
       _E_INDU      = E_U;
   }

   void setEnergy_PPUU(double epp, 
                       double epu, 
                       double euu) {

       _EPP = epp;
       _EPU = epu;
       _EUU = euu;           
   }

   void  setEnergy_f_m(double e_f_c_c,
                       double e_f_c_non_c,
                       double e_f_c_out,
                       double e_f_non_c_non_c, 
                       double e_f_non_c_out,
                       double e_m_c,       
                       double e_m_non_c,                              
                       double e_m_out) {

       _EF_PAIR_PAIR   = e_f_c_c;
       _EF_PAIR_SPH1   = e_f_c_non_c;
       _EF_PAIR_SPH2   = e_f_c_out;
       _EF_SPH1_SPH1   = e_f_non_c_non_c;
       _EF_SPH1_SPH2   = e_f_non_c_out;

       _EM_PAIR        = e_m_c;
       _EM_SPH1        = e_m_non_c;
       _EM_SPH2        = e_m_out;
   }
   
   void setEnergy_QMMM(double energy_QM,
                       double energy_GWBSE,
                       double energy_SF, 
                       double energy_QMMM) {
       
       _E_QM = energy_QM;
       _E_SF = energy_SF;
       _E_GWBSE = energy_GWBSE;
       _E_QMMM = energy_QMMM;
   }

private:
  
   int                  _id;
   int                  _userId;
   std::string               _tag;
   Topology            *_top;   
   
   bool                 _start_from_cpt;
   std::vector<Segment*>     _qmSegs;
   std::vector<std::string>       _qmSegMps;
   vec                  _center;
   std::map<int,bool>        _isSegInCenter;
   PolarTop            *_ptop;
   bool                 _clean_ptop;

   // Inductor facts
   int                  _qm0_size;
   int                  _mm1_size;
   int                  _mm2_size;   
   int                  _iter;

   // Energy splittings:
   double               _E_Tot;
   // ... 0th kind
   double               _E_Pair_Pair;
   double               _E_Pair_Sph1;
   double               _E_Pair_Sph2;
   double               _E_Sph1_Sph1;
   double               _E_Sph1_Sph2;
   // ... 1st kind
   double               _E_PERM;
   double               _E_INDU;
   // ... 2nd kind
   double               _EPP;
   double               _EPU;
   double               _EUU;
   // ... 3rd kind
   double               _EF_PAIR_PAIR;
   double               _EF_PAIR_SPH1;
   double               _EF_PAIR_SPH2;
   double               _EF_SPH1_SPH1;
   double               _EF_SPH1_SPH2;
   double               _EM_PAIR;
   double               _EM_SPH1;
   double               _EM_SPH2;  
   
   // QM Energy (self-energy, internal)
   double               _E_QM;
   double               _E_SF;
   double               _E_GWBSE;
   double               _E_QMMM;
   
   std::string               _infoLine;

}; 
    

inline bool XJob::isInCenter(int segId) {

   bool inCenter = (this->_isSegInCenter.count(segId) > 0) ? true : false;
   return inCenter;
}

inline bool XJob::isWithinDist(const vec &pt, double dist, Topology *top) {
//    // This will not work for pairs far apart
//    bool yesno = false;
//    
//    double dR = abs(top->PbShortestConnect(_center, pt));
//    if (dR <= dist) { yesno = true; }
//    
//    return yesno;
    
    bool inCenter = false;
    
    for (unsigned int i = 0; i < _qmSegs.size(); ++i) {         
        Segment *seg = _qmSegs[i];
        double dR = abs(top->PbShortestConnect(seg->getPos(), pt));
        if (dR <= dist) { inCenter = true; break; }
    }
    
    return inCenter;    
}

template<typename JobContainer, typename pJob>
JobContainer XJOBS_FROM_TABLE(const std::string &job_file, Topology *top);
    
}}

#endif // VOTCA_XTP_XJOB_H

