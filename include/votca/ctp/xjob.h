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

namespace votca { namespace ctp {
  
// +++++++++++++++++++++++++++ //
// XJob (distr. among threads) //
// +++++++++++++++++++++++++++ //

class XJob
{
public:

    XJob(int id,      string tag,  string job_type, int site_id,
         int pairId,  int seg1Id,  int seg2Id,
         string mps1, string mps2, Topology *top) :

         _id(id),         _tag(tag),       _pairId(pairId),
         _type(job_type), _seg1Id(seg1Id), _seg2Id(seg2Id),
         _mps1(mps1),     _mps2(mps2),     _start_from_cpt(false),
         _site_id(site_id) {

         _seg1 = top->getSegment(seg1Id);
         _seg2 = top->getSegment(seg2Id);      

         // JOB TYPE: PAIR
         if (job_type == "pair") {

             // Watch out for periodic-boundary correction
             _pair = top->NBList().FindPair(_seg1,_seg2);

             if (_pair == NULL) {
                cout << endl
                     << "WARNING: No such pair " << pairId
                     << "(SEG1 " << seg1Id << ", SEG2 " << seg2Id << "). "
                     << flush;
                throw runtime_error("Pair specs do not match topology.");
             }

             assert ( pairId == _pair->getId() );                           // Inconsistency in job input w.r.t. top

             _center = ( _pair->Seg1PbCopy()->Atoms().size()
                       * _pair->Seg1PbCopy()->getPos()

                       + _pair->Seg2PbCopy()->Atoms().size()
                       * _pair->Seg2PbCopy()->getPos() )

                       / (_pair->Seg1PbCopy()->Atoms().size()
                       +  _pair->Seg2PbCopy()->Atoms().size() );
         }

         // JOB TYPE: SITE
         else if (job_type == "site") {

             if      (_site_id == seg1Id) { _center = _seg1->getPos();
                                            _site   = _seg1;           }
             else if (_site_id == seg2Id) { _center = _seg2->getPos();
                                            _site   = _seg2;           }
             else    { throw runtime_error("ERROR: "
                   "No such site " + boost::lexical_cast<string>(_site_id) +
                   " in job " + boost::lexical_cast<string>(_id) + ". "); }
         }

         // JOB TYPE: UNKNOWN
         else {
            cout << endl 
                 << "ERROR: Job " << _id << ": Invalid job type "
                 << job_type << ": Should be either 'pair' or 'site'"
                 << endl;
            throw runtime_error("ERROR: Faulty input in job list.");
         }
   }

   ~XJob() {};

   int      getId()                     { return _id; }
   string   getTag()                    { return _tag; }
   int      getPairId()                 { return _pairId; }
   int      getSeg1Id()                 { return _seg1Id; }
   int      getSeg2Id()                 { return _seg2Id; }
   string   getMPS1()                   { return _mps1; }
   string   getMPS2()                   { return _mps2; }
   int      getSiteId()                 { return _site_id; }
   string   getType()                   { return _type; }

   void     setType(string typ, int id) { _type = typ; _site_id = id; }
   void     setIter(int iter)           { _iter = iter; }
   void     setSizePol(int size)        { _sizePol = size; }
   void     setSizeShell(int size)      { _sizeShell = size; }

   bool     isInCenter(int segId) {

       bool inCenter = false;

       // Job type 'pair'
       if (this->_type == "pair") {               
           if (segId == this->_seg1Id || segId == this->_seg2Id) {
               inCenter = true;
           }
           else {
               inCenter = false;
           }
       }
       // Job type 'site'
       else if (this->_type == "site") {
           if (segId == this->_site_id) {
               inCenter = true;
           }
           else {
               inCenter = false;
           }
       }
       // Unrecognised job type
       else {
           assert(false);
       }

       return inCenter;
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

   void     WriteInfoLine(FILE *out) {

       fprintf(out, "%5d %-20s  ",
                    _id, _tag.c_str() );
       fprintf(out, "E_TOT %+4.7f ",
                    _E_Tot );
       fprintf(out, "EPP %+4.7f "
                    "EPU %+4.7f "
                    "EUU %+4.7f ",
                    _EPP,
                    _EPU,
                    _EUU);
       fprintf(out, "EF0_0 %+4.7f "
                    "EF0_1 %+4.7f "
                    "EF0_2 %+4.7f "
                    "EF1_1 %+4.7f "
                    "EF1_2 %+4.7f "
                    "EM_0_ %+4.7f "
                    "EM_1_ %+4.7f "
                    "EM_2_ %+4.7f ",
                    _EF_PAIR_PAIR,
                    _EF_PAIR_SPH1,
                    _EF_PAIR_SPH2,
                    _EF_SPH1_SPH1, 
                    _EF_SPH1_SPH2,
                    _EM_PAIR,
                    _EM_SPH1,      
                    _EM_SPH2);
       fprintf(out, "ITER %3d SPHERE %4d SHELL %4d "
                    "CENTER %4.7f %4.7f %4.7f ",
                    _iter, _sizePol, _sizeShell, 
                    _center.getX(), _center.getY(), _center.getZ() );
       fprintf(out, "E_CUT0_CUT0 %+4.7f "
                    "E_CUT0_CUT1 %+4.7f "
                    "E_CUT0_CUT2 %+4.7f "
                    "E_CUT1_CUT1 %+4.7f "
                    "E_CUT1_CUT2 %+4.7f "
                    "E_PERM %+4.7f "
                    "E_INDU %+4.7f ",
                    _E_Pair_Pair,
                    _E_Pair_Sph1,
                    _E_Pair_Sph2,
                    _E_Sph1_Sph1,  
                    _E_Sph1_Sph2,
                    _E_PERM,
                    _E_INDU );
       fprintf(out, "\n");

   }

   bool     StartFromCPT()  { return _start_from_cpt; }

   Segment *Seg1()          { return _seg1; }
   Segment *Seg2()          { return _seg2; }
   QMPair  *Pair()          { return _pair; }
   vec     &Center()        { return _center; }

private:

   int          _id;
   string       _tag;
   string       _type;
   int          _site_id; // only relevant if type == "site"
   int          _pairId;
   int          _seg1Id;
   int          _seg2Id;
   string       _mps1;
   string       _mps2;

   bool         _start_from_cpt;

   QMPair      *_pair;
   Segment     *_site;
   Segment     *_seg1;
   Segment     *_seg2;
   vec          _center;


   int          _iter;
   int          _sizePol;
   int          _sizeShell;

   // Energy splittings:
   double       _E_Tot;
   // ... 0th kind
   double       _E_Pair_Pair;
   double       _E_Pair_Sph1;
   double       _E_Pair_Sph2;
   double       _E_Sph1_Sph1;
   double       _E_Sph1_Sph2;
   // ... 1st kind
   double       _E_PERM;
   double       _E_INDU;
   // ... 2nd kind
   double       _EPP;
   double       _EPU;
   double       _EUU;
   // ... 3rd kind
   double       _EF_PAIR_PAIR;
   double       _EF_PAIR_SPH1;
   double       _EF_PAIR_SPH2;
   double       _EF_SPH1_SPH1;
   double       _EF_SPH1_SPH2;
   double       _EM_PAIR;
   double       _EM_SPH1;
   double       _EM_SPH2;       

}; 
    

vector<XJob*> XJOBS_FROM_TABLE(const string &job_file, Topology *top);


    
}}

#endif