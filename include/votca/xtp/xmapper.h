/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __XMAPPER__H
#define	__XMAPPER__H

#include <votca/tools/mutex.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/apolarsite.h>
#include <votca/xtp/qmthread.h>

// TODO Change maps to _alloc_xmlfile_fragsegmol_***
// TODO Confirm thread safety
// TODO Add "const" keyword

namespace votca { namespace xtp {
    
class XMpsMap
{

public:        

    XMpsMap() : _alloc_table("no_alloc"), _estatics_only(false) {};
   ~XMpsMap() {};

    // User interface:
    void GenerateMap(string xml_file, string alloc_table, Topology *top);
    void EquipWithPolSites(Topology *top);
    
    // Adapt to XJob
    PolarSeg *MapPolSitesToSeg(const vector<APolarSite*> &pols_n, Segment *seg, bool only_active_sites = true);
    vector<APolarSite*> GetOrCreateRawSites(const string &mpsfile, QMThread *thread = NULL);
    void Gen_QM_MM1_MM2(Topology *top, XJob *job, double co1, double co2, QMThread *thread = NULL);
    void Gen_FGC_FGN_BGN(Topology *top, XJob *job, QMThread *thread = NULL);
    void Gen_BGN(Topology *top, PolarTop *ptop, QMThread *thread = NULL);
    void Gen_FGC_Load_FGN_BGN(Topology *top, XJob *job, string archfile, QMThread *thread = NULL);
    
    void setEstaticsOnly(bool estatics_only) { _estatics_only = estatics_only; }
    
    // Called by GenerateMap(...)
    void CollectMapFromXML(string xml_file);
    void CollectSegMpsAlloc(string alloc_table, Topology *top);
    void CollectSitesFromMps();
    
    
private:

    string _alloc_table;
    votca::tools::Mutex  _lockThread;
    bool _estatics_only;
    
    // Maps retrieved from XML mapping files
    map<string, bool>                   _map2md;
    map<string, vector<int> >           _alloc_frag_mpoleIdx;
    map<string, vector<string> >        _alloc_frag_mpoleName;
    map<string, vector<int> >           _alloc_frag_trihedron;
    map<string, vector<double> >        _alloc_frag_weights;
    map<string, vector<int> >           _alloc_frag_isVirtualMp;

    // Allocation of mps-files to segments, state-resolved
    map<int,string>                 _segId_mpsFile_n;
    map<int,string>                 _segId_mpsFile_e;
    map<int,string>                 _segId_mpsFile_h;

    // Raw polar sites collected from mps-files
    map<string,vector<APolarSite*> > _mpsFile_pSites;
    map<string,vector<APolarSite*> > _mpsFile_pSites_job;
};
    
    

    
}}


#endif