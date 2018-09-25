/*
 *            Copyright 2016 The MUSCET Development Team
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


#include <votca/xtp/qmmolecule.h>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>
#include "qmmolecule.h"


using namespace votca::tools;

namespace votca { namespace xtp {

    QMMolecule::ReadFromCpt(CptLoc parent) {

        CptLoc qmAtomsGr = parent.openGroup("qmatoms");
        size_t count = qmAtomsGr.getNumObjs();
      _atomlist.resize(0);
      _atomlist.reserve(count);
        for (size_t i = 0; i < count; ++i) {
        CptLoc tempLoc = qmAtomsGr.openGroup("atom" + std::to_string(i));
        QMAtom temp;
        temp.ReadFromCpt(tempLoc);
        _atomlist.push_back(temp);
      }
      

    }

    QMMolecule::WriteToCpt(CptLoc parent) const {
      CptLoc qmAtomsGr = parent.createGroup("qmatoms");
      size_t count = 0;
      for (const auto& qma : _atomlist) {
        CptLoc tempLoc = qmAtomsGr.createGroup("atom" + std::to_string(count));
        qma->WriteToCpt(tempLoc);
        ++count;
      }

    }
    

}}


