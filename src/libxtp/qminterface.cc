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
                                                                                           
#include <votca/xtp/qmatom.h>                                                                                                      
#include <votca/xtp/qminterface.h>
#include <votca/xtp/segment.h>  
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/atom.h>
#include <boost/format.hpp>

using boost::format;

namespace votca {
    namespace xtp {

        QMMolecule QMInterface::Convert(std::vector<Segment* > segments) {
            QMMolecule result=QMMolecule("QMAtoms",0);
            for (unsigned AtomId=0;AtomId<segments.size();++AtomId) {
                Segment* segment=segments[AtomId];
                std::vector < Atom* >& atoms = segment->Atoms();
                for (Atom* atom: atoms) {
                    if(!atom->HasQMPart()){
                      continue;
                    }
                    Eigen::Vector3d pos = atom->getQMPos().toEigen() * tools::conv::nm2bohr;
                    result.push_back(QMAtom(AtomId,atom->getElement(),pos));
                }
            }
            return result;
        }

       

   

    }
}
