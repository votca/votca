/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_MAPCHECKER_H
#define _VOTCA_XTP_MAPCHECKER_H


#include <votca/xtp/qmcalculator.h>


namespace votca {
namespace xtp {

class MapChecker : public QMCalculator {
 public:
  MapChecker(){};

  ~MapChecker(){};

  std::string Identify() { return "mapchecker"; }

  void Initialize(tools::Property &opt);
  bool EvaluateFrame(Topology &top);
 

 private:
     std::string _segmentfile;
     std::string _qmfile;
     std::string _mpfile;

};

void MapChecker::Initialize(tools::Property &opt) {
    std::string key = "options." + Identify();
_segmentfile = opt.ifExistsReturnElseReturnDefault<std::string>(
      key + ".segment_pdbfile", "segments.pdb");

_qmfile = opt.ifExistsReturnElseReturnDefault<std::string>(
      key + ".qm_pdbfile", "qm_segments.pdb");

_mpfile = opt.ifExistsReturnElseReturnDefault<std::string>(
      key + ".mp_pdbfile", "mp_segments.pdb");

  
}

bool MapChecker::EvaluateFrame(Topology &top) {

   std::string base=tools::filesystem::GetFileBase(_segmentfile);
   std::string fileending=tools::filesystem::GetFileExtension(_segmentfile);
   top.WriteToPdb(base+"_"+top.getStep()+"."+fileending);

  return true;
}


}  // namespace xtp
}  // namespace votca

#endif  // _VOTCA_XTP_MAPCHECKER_H
