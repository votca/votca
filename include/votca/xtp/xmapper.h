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

#ifndef VOTCA_XTP_XMAPPER_H
#define	VOTCA_XTP_XMAPPER_H

#include <votca/tools/mutex.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/qmthread.h>

// TODO Change maps to _alloc_xmlfile_fragsegmol_***
// TODO Confirm thread safety
// TODO Add "const" keyword

namespace votca { namespace xtp {

class Topology;
class Segment;
class APolarSite;
class PolarSeg;
 
class XMpsMap
{

public:        

    XMpsMap() : _alloc_table("no_alloc"), _estatics_only(false) {};
   ~XMpsMap() {};

    // User interface:
    void GenerateMap(std::string xml_file, std::string alloc_table, Topology *top);
    void EquipWithPolSites(Topology *top);
    
    // Adapt to XJob
    PolarSeg *MapPolSitesToSeg(const std::vector<APolarSite*> &pols_n, Segment *seg, bool only_active_sites = true);
    std::vector<APolarSite*> GetOrCreateRawSites(const std::string &mpsfile, QMThread *thread = NULL);
    void Gen_QM_MM1_MM2(Topology *top, XJob *job, double co1, double co2, QMThread *thread = NULL);
    void Gen_FGC_FGN_BGN(Topology *top, XJob *job, QMThread *thread = NULL);
    void Gen_BGN(Topology *top, PolarTop *ptop, QMThread *thread = NULL);
    void Gen_FGC_Load_FGN_BGN(Topology *top, XJob *job, std::string archfile, QMThread *thread = NULL);
    
    void setEstaticsOnly(bool estatics_only) { _estatics_only = estatics_only; }
    
    // Called by GenerateMap(...)
    void CollectMapFromXML(std::string xml_file);
    void CollectSegMpsAlloc(std::string alloc_table, Topology *top);
    void CollectSitesFromMps();
    
    
private:

    std::string _alloc_table;
    votca::tools::Mutex  _lockThread;
    bool _estatics_only;
    
    // Maps retrieved from XML mapping files
    std::map<std::string, bool>                   _map2md;
    std::map<std::string, std::vector<int> >           _alloc_frag_mpoleIdx;
    std::map<std::string, std::vector<std::string> >        _alloc_frag_mpoleName;
    std::map<std::string, std::vector<int> >           _alloc_frag_trihedron;
    std::map<std::string, std::vector<double> >        _alloc_frag_weights;
    std::map<std::string, std::vector<int> >           _alloc_frag_isVirtualMp;

    // Allocation of mps-files to segments, state-resolved
    std::map<int,std::string>                 _segId_mpsFile_n;
    std::map<int,std::string>                 _segId_mpsFile_e;
    std::map<int,std::string>                 _segId_mpsFile_h;

    // Raw polar sites collected from mps-files
    std::map<std::string,std::vector<APolarSite*> > _mpsFile_pSites;
    std::map<std::string,std::vector<APolarSite*> > _mpsFile_pSites_job;
};
    
    

    
}}


#endif // VOTCA_XTP_XMAPPER_H

