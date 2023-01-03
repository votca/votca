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

#ifndef __XMAPPER__H
#define __XMAPPER__H

#include "votca/xtp/backgroundregion.h"
#include "votca/xtp/segmentmapper.h"
#include <votca/tools/mutex.h>
#include <votca/xtp/ewald/apolarsite.h>
#include <votca/xtp/ewald/xjob.h>
#include <votca/xtp/topology.h>

// TODO Change maps to _alloc_xmlfile_fragsegmol_***
// TODO Confirm thread safety
// TODO Add "const" keyword

namespace votca {
namespace xtp {

class XMapper {

 public:
  // Adapt to XJob
  void setLogger(Logger* log) { log_ = log; };
  // void Gen_QM_MM1_MM2(Topology *top, XJob *job, double co1, double co2,
  // QMThread *thread = NULL);
  void Gen_FGC_FGN_BGN(std::string mapfile, const Topology& top, XJob* xjob);
  // void Gen_BGN(Topology *top, PolarTop *ptop, QMThread *thread = NULL);
  void Gen_FGC_Load_FGN_BGN(std::string mapfile, const Topology& top, XJob* job,
                            std::string archfile);
  // QMThread *thread = NULL);

 private:
  // logger
  Logger* log_;

  std::string _alloc_table;
  votca::tools::Mutex _lockThread;

  // Maps retrieved from XML mapping files
  std::map<std::string, bool> _map2md;
  std::map<std::string, std::vector<int> > _alloc_frag_mpoleIdx;
  std::map<std::string, std::vector<std::string> > _alloc_frag_mpoleName;
  std::map<std::string, std::vector<int> > _alloc_frag_trihedron;
  std::map<std::string, std::vector<double> > _alloc_frag_weights;
  std::map<std::string, std::vector<int> > _alloc_frag_isVirtualMp;

  // Allocation of mps-files to segments, state-resolved
  std::map<int, std::string> _segId_mpsFile_n;
  std::map<int, std::string> _segId_mpsFile_e;
  std::map<int, std::string> _segId_mpsFile_h;

  // Raw polar sites collected from mps-files
  std::map<std::string, std::vector<APolarSite*> > _mpsFile_pSites;
  std::map<std::string, std::vector<APolarSite*> > _mpsFile_pSites_job;
};

}  // namespace xtp
}  // namespace votca

#endif