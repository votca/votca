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




#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/elements.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/qminterface.h>

#include "votca/xtp/qmatom.h"

using boost::format;

namespace votca {
    namespace xtp {

        ctp::APolarSite *QMInterface::Convert(QMAtom *atm, int id) {
          
            
            std::string elem = atm->getType();
            ctp::APolarSite *new_aps = new ctp::APolarSite(id, elem);
            
            double pol = 0.0;
            try {
                pol = _polar_table.at(elem);
            } catch (const std::exception& out_of_range) {
                std::cout << std::endl << "QMMInterface - no default polarizability given "
                        << "for element type '" << elem << "'. Defaulting to 1A**3" << std::flush;
                pol = 1e-3;
            }
            new_aps->setIsoP(pol);
                   
            tools::vec pos = atm->getPos()*tools::conv::bohr2nm;
            double q = atm->getPartialcharge();
            new_aps->setRank(0);
            new_aps->setPos(pos);
            new_aps->setQ00(q, 0); // <- charge state 0 <> 'neutral'
           
            return new_aps;
        }

        ctp::PolarSeg QMInterface::Convert(std::vector<QMAtom*> &atms) {
            ctp::PolarSeg new_pseg = ctp::PolarSeg();
            std::vector<QMAtom*>::iterator it;
            for (it = atms.begin(); it < atms.end(); ++it) {
                ctp::APolarSite *new_site = this->Convert(*it);
                new_pseg.push_back(new_site);
            }
            return new_pseg;
        }

        std::vector<QMAtom *> QMInterface::Convert(std::vector<ctp::Segment* > segments) {

            std::vector<QMAtom *> qmatoms;

            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
            int AtomId=0;
            for (sit = segments.begin(); sit != segments.end(); ++sit) {
                std::vector < ctp::Atom* >& _atoms = (*sit)->Atoms();

                for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                    if(!(*ait)->HasQMPart()){
                      continue;
                    }

                    tools::vec pos = (*ait)->getQMPos() * tools::conv::nm2bohr;
                    std::string name = (*ait)->getElement();
                    //be careful charges are set to zero because most of the time ctp::Atom has getQ not set, the construct is very weird, ctp is shit
                    QMAtom* qmatom = new QMAtom(AtomId,name, pos);
                    AtomId++;
                    qmatoms.push_back(qmatom);

                }
            }
            return qmatoms;
        }

        void QMInterface::GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop, Orbitals &orb) {
          int AtomID=0;
            // INNER SHELL QM0
            for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->QM0()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    ctp::APolarSite *aps = (*pseg)[j];
                    tools::vec pos = aps->getPos() * tools::conv::nm2bohr;
                    orb.AddAtom(AtomID,aps->getName(), pos);
                    AtomID++;
                }
            }

            return;
        }


        std::vector<ctp::PolarSeg*> QMInterface::GenerateMultipoleList(ctp::PolarTop *ptop  ) {

            std::vector<ctp::PolarSeg*> MultipoleList;

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                ctp::PolarSeg *pseg = new ctp::PolarSeg(ptop->MM1()[i],false);
                MultipoleList.push_back(pseg);
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                ctp::PolarSeg *pseg = new ctp::PolarSeg(ptop->MM2()[i],false);
                MultipoleList.push_back(pseg);
            }

            return MultipoleList;
        }


        void QMInterface::Orbitals2Segment(ctp::Segment* _segment, const Orbitals& _orbitals) {

            
            std::vector<QMAtom* > ::iterator ait;
            std::vector< QMAtom* >_atoms = _orbitals.QMAtoms();

            std::string type;
            int id = 1;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {


                type = (*ait)->getType();
                
                ctp::Atom *pAtom = new ctp::Atom(id++, type);
                tools::vec position= (*ait)->getPos()*votca::tools::conv::bohr2nm;
                pAtom->setPos(position);
                pAtom->setQMPart(id, position);
                pAtom->setElement(type);
                _segment->AddAtom(pAtom);
            }
            return;
        }


    }
}
