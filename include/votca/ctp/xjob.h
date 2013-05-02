/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __XJOB__H
#define	__XJOB__H

#include <votca/ctp/topology.h>
#include <votca/ctp/polartop.h>

namespace votca { namespace ctp {
  

class XJob
{
public:

    XJob(int id, string tag, vector<Segment*> &qmSegs, 
         vector<string> &qmSegMps, Topology *top);

   ~XJob();

   int                  getId()                     { return _id; }
   string               getTag()                    { return _tag; }
   vector<Segment*>    &getSegments()               { return _qmSegs; }
   vector<string>      &getSegMps()                 { return _qmSegMps; }
   Topology            *getTop()                    { return _top; }
   PolarTop            *getPolarTop()               { return _ptop; }
   int                  getIter()                   { return _iter; }   
   void                 setPolarTop(PolarTop *ptop);
   
   void                 CalcCenterPos();
   vec                 &Center()                    { return _center; }
   
   inline bool          isInCenter(int segId);
   inline bool          isWithinDist(const vec &pt, double dist, Topology *top);
   bool                 StartFromCPT()  { return _start_from_cpt; }
   void                 WriteInfoLine(FILE *out);

   
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

private:
  
   int                  _id;
   string               _tag;
   Topology            *_top;   
   
   vector<string>       _qmSegMps;
   vector<Segment*>     _qmSegs;
   vec                  _center;
   map<int,bool>        _isSegInCenter;
   bool                 _start_from_cpt;
   PolarTop            *_ptop;

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

}; 
    

inline bool XJob::isInCenter(int segId) {

   bool inCenter = (this->_isSegInCenter.count(segId) > 0) ? true : false;
   return inCenter;
}

inline bool XJob::isWithinDist(const vec &pt, double dist, Topology *top) {
    bool yesno = false;
    
    double dR = abs(top->PbShortestConnect(_center, pt));
    if (dR <= dist) { yesno = true; }
    
    return yesno;
}


vector<XJob*> XJOBS_FROM_TABLE(const string &job_file, Topology *top);


    
}}

#endif