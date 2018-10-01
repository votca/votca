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

        APolarSite *QMInterface::Convert(QMAtom *atm, int id) {
                
            std::string elem = atm->getElement();
            APolarSite *new_aps = new APolarSite(id, elem);
            
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

        PolarSeg QMInterface::Convert(std::vector<QMAtom*> &atms) {
            PolarSeg new_pseg = PolarSeg();
            std::vector<QMAtom*>::iterator it;
            for (it = atms.begin(); it < atms.end(); ++it) {
                APolarSite *new_site = this->Convert(*it);
                new_pseg.push_back(new_site);
            }
            return new_pseg;
        }

        std::vector<QMAtom *> QMInterface::Convert(std::vector<Segment* > segments) {
            std::vector<QMAtom *> qmatoms;
            int AtomId=0;
            for (Segment* segment:segments) {
                std::vector < Atom* >& atoms = segment->Atoms();
                for (Atom* atom: atoms) {
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

        QMMolecule QMInterface::GenerateQMAtomsFromPolarSegs(PolarTop *ptop) {
          int AtomID=0;
            for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
                PolarSeg *pseg = ptop->QM0()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    APolarSite *aps = (*pseg)[j];
                    tools::vec pos = aps->getPos() * tools::conv::nm2bohr;
                    orb.AddAtom(AtomID,aps->getName(), pos);
                    AtomID++;
                }
            }

            return;
        }


        std::vector<std::shared_ptr<PolarSeg> > QMInterface::GenerateMultipoleList(PolarTop *ptop  ) {
            std::vector<std::shared_ptr<PolarSeg> > MultipoleList;

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                std::shared_ptr<PolarSeg>  pseg ( new PolarSeg(ptop->MM1()[i],false));
                MultipoleList.push_back(pseg);
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                std::shared_ptr<PolarSeg>  pseg (new PolarSeg(ptop->MM2()[i],false));
                MultipoleList.push_back(pseg);
            }
            return MultipoleList;
        }

        void QMInterface::Orbitals2Segment(Segment& segment, const Orbitals& orbitals) {

            std::vector< QMAtom* >_atoms = orbitals.QMAtoms();
            std::string type;
            int id = 1;
            for (QMAtom* atom: _atoms) {
                type = atom->getType();
                Atom *pAtom = new Atom(id++, type);
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
