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
#include "votca/xtp/atom.h"

using boost::format;

namespace votca {
    namespace xtp {

        xtp::APolarSite *QMMInterface::Convert(QMAtom *atm, int id) {
          
            
            std::string elem = atm->getType();
            xtp::APolarSite *new_aps = new xtp::APolarSite(id, elem);
            
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

        xtp::PolarSeg QMMInterface::Convert(std::vector<QMAtom*> &atms) {
            xtp::PolarSeg new_pseg = xtp::PolarSeg();
            std::vector<QMAtom*>::iterator it;
            for (it = atms.begin(); it < atms.end(); ++it) {
                xtp::APolarSite *new_site = this->Convert(*it);
                new_pseg.push_back(new_site);
            }
            return new_pseg;
        }

        std::vector<QMAtom *> QMMInterface::Convert(std::vector<xtp::Segment* > segments) {

            std::vector<QMAtom *> qmatoms;

            std::vector< xtp::Atom* > ::iterator ait;
            std::vector< xtp::Segment* >::iterator sit;
            int AtomId=0;
            for (sit = segments.begin(); sit != segments.end(); ++sit) {
                std::vector < xtp::Atom* >& _atoms = (*sit)->Atoms();

                for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                    if(!(*ait)->HasQMPart()){
                      continue;
                    }

                    tools::vec pos = (*ait)->getQMPos() * tools::conv::nm2bohr;
                    std::string name = (*ait)->getElement();
                    //be careful charges are set to zero because most of the time xtp::Atom has getQ not set, the construct is very weird, xtp is shit
                    QMAtom* qmatom = new QMAtom(AtomId,name, pos);
                    AtomId++;
                    qmatoms.push_back(qmatom);

                }
            }
            return qmatoms;
        }

        void QMMInterface::GenerateQMAtomsFromPolarSegs(xtp::PolarTop *ptop, Orbitals &orb) {
          int AtomID=0;
            // INNER SHELL QM0
            for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
                xtp::PolarSeg *pseg = ptop->QM0()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    xtp::APolarSite *aps = (*pseg)[j];
                    tools::vec pos = aps->getPos() * tools::conv::nm2bohr;
                    orb.AddAtom(AtomID,aps->getName(), pos);
                    AtomID++;
                }
            }

            return;
        }


        std::vector<xtp::PolarSeg*> QMMInterface::GenerateMultipoleList(xtp::PolarTop *ptop  ) {

            std::vector<xtp::PolarSeg*> MultipoleList;

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                xtp::PolarSeg *pseg = new xtp::PolarSeg(ptop->MM1()[i],false);
                MultipoleList.push_back(pseg);
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                xtp::PolarSeg *pseg = new xtp::PolarSeg(ptop->MM2()[i],false);
                MultipoleList.push_back(pseg);
            }

            return MultipoleList;
        }


        void QMMInterface::Orbitals2Segment(xtp::Segment* _segment, Orbitals* _orbitals) {

            
            std::vector<QMAtom* > ::iterator ait;
            std::vector< QMAtom* >_atoms = _orbitals->QMAtoms();

            std::string type;
            int id = 1;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {


                type = (*ait)->getType();
                
                xtp::Atom *pAtom = new xtp::Atom(id++, type);
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
