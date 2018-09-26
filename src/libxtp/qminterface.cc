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

#include <votca/tools/elements.h>                                                

#include <votca/xtp/apolarsite.h>                                                
#include <votca/xtp/atom.h>
#include <votca/xtp/qmatom.h>                                                    
#include <votca/xtp/orbitals.h>   
#include <votca/xtp/polarseg.h>                                                  
#include <votca/xtp/polartop.h>                                                  
#include <votca/xtp/qminterface.h>
#include <votca/xtp/segment.h>                                                   

using boost::format;

namespace votca {
    namespace xtp {

        xtp::APolarSite *QMInterface::Convert(QMAtom *atm, int id) {
                
            std::string elem = atm->getElement();
            xtp::APolarSite *new_aps = new xtp::APolarSite(id, elem);
            
            double pol = 0.0;
            try {
                pol = _element.getPolarizability(elem);
            } catch (const std::invalid_argument& out_of_range) {
                std::cout << std::endl << "QMInterface - no default polarizability given "
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

        xtp::PolarSeg QMInterface::Convert(std::vector<QMAtom*> &atms) {
            xtp::PolarSeg new_pseg = xtp::PolarSeg();
            std::vector<QMAtom*>::iterator it;
            for (it = atms.begin(); it < atms.end(); ++it) {
                xtp::APolarSite *new_site = this->Convert(*it);
                new_pseg.push_back(new_site);
            }
            return new_pseg;
        }

        std::vector<QMAtom *> QMInterface::Convert(std::vector<xtp::Segment* > segments) {
            std::vector<QMAtom *> qmatoms;
            int AtomId=0;
            for (xtp::Segment* segment:segments) {
                std::vector < xtp::Atom* >& atoms = segment->Atoms();
                for (xtp::Atom* atom: atoms) {
                    if(!atom->HasQMPart()){
                      continue;
                    }
                    tools::vec pos = atom->getQMPos() * tools::conv::nm2bohr;
                    std::string name = atom->getElement();
                    QMAtom* qmatom = new QMAtom(AtomId,name, pos);
                    AtomId++;
                    qmatoms.push_back(qmatom);
                }
            }
            return qmatoms;
        }

        void QMInterface::GenerateQMAtomsFromPolarSegs(xtp::PolarTop *ptop, Orbitals &orb) {
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


        std::vector<std::shared_ptr<xtp::PolarSeg> > QMInterface::GenerateMultipoleList(xtp::PolarTop *ptop  ) {
            std::vector<std::shared_ptr<xtp::PolarSeg> > MultipoleList;

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                std::shared_ptr<xtp::PolarSeg>  pseg ( new xtp::PolarSeg(ptop->MM1()[i],false));
                MultipoleList.push_back(pseg);
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                std::shared_ptr<xtp::PolarSeg>  pseg (new xtp::PolarSeg(ptop->MM2()[i],false));
                MultipoleList.push_back(pseg);
            }
            return MultipoleList;
        }

        void QMInterface::Orbitals2Segment(xtp::Segment& segment, const Orbitals& orbitals) {

            std::vector< QMAtom* >_atoms = orbitals.QMAtoms();
            std::string type;
            int id = 1;
            for (QMAtom* atom: _atoms) {
                type = atom->getType();
                xtp::Atom *pAtom = new xtp::Atom(id++, type);
                tools::vec position= atom->getPos()*votca::tools::conv::bohr2nm;
                pAtom->setPos(position);
                pAtom->setQMPart(id, position);
                pAtom->setElement(type);
                segment.AddAtom(pAtom);
            }
            return;
        }


    }
}
