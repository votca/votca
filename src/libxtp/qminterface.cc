/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>


#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/xtp/elements.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/qminterface.h>

using boost::format;

namespace votca {
    namespace xtp {

    ctp::APolarSite *QMMInterface::Convert(ctp::QMAtom *atm, int id) {
        double A_to_nm = 0.1;
        vec pos = A_to_nm*vec(atm->x, atm->y, atm->z);
        double q = atm->charge;
        std::string elem = atm->type;
        double pol = 0.0;
        try {
            pol = _polar_table.at(elem);
        }
        catch(out_of_range) {
            std::cout << std::endl << "QMMInterface - no default polarizability given "
                << "for element type '" << elem << "'. Defaulting to 1A**3" << std::flush;
            pol = 1e-3;
        }

        ctp::APolarSite *new_aps = new ctp::APolarSite(id, elem);
        new_aps->setRank(0);
        new_aps->setPos(pos);
        new_aps->setQ00(q,0); // <- charge state 0 <> 'neutral'
        new_aps->setIsoP(pol);
        
        return new_aps;
    }
        
    ctp::PolarSeg QMMInterface::Convert(std::vector<ctp::QMAtom*> &atms) {        
        ctp::PolarSeg new_pseg =ctp::PolarSeg();
        std::vector<ctp::QMAtom*>::iterator it;
        for (it = atms.begin(); it < atms.end(); ++it) {
            ctp::APolarSite *new_site = this->Convert(*it);
            new_pseg.push_back(new_site);
        }
        return new_pseg;
    }
    
    
    std::vector<ctp::QMAtom *> QMMInterface::Convert( std::vector<ctp::Segment* > segments){
        
        std::vector<ctp::QMAtom *> qmatoms;
        
        std::vector< ctp::Atom* > ::iterator ait;
        std::vector< ctp::Segment* >::iterator sit;
        for (sit = segments.begin(); sit != segments.end(); ++sit) {
            std::vector < ctp::Atom* >& _atoms = (*sit)->Atoms();

            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {


                vec pos = (*ait)->getQMPos()*tools::conv::nm2ang;
                std::string name = (*ait)->getElement();
                //be careful charges are set to zero because most of the time ctp::Atom has getQ not set, the construct is very weird, ctp is shit
                ctp::QMAtom* qmatom=new ctp::QMAtom(name,pos,0.0,!(*ait)->HasQMPart());
               
                qmatoms.push_back(qmatom);
                        
            }
        }
        return qmatoms;
    }
    
    
    void QMMInterface::GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop, Orbitals &orb) {

            // INNER SHELL QM0
            for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->QM0()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    ctp::APolarSite *aps = (*pseg)[j];
                    vec pos = aps->getPos() *tools::conv::nm2ang;
                    double Q = aps->getQ00();
                    orb.AddAtom(aps->getName(), pos, Q, false);
                }
            }

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->MM1()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    ctp::APolarSite *aps = (*pseg)[j];
                    addMMAtomtoOrb(aps,orb,true);
                }
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->MM2()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {
                    ctp::APolarSite *aps = (*pseg)[j];
                    addMMAtomtoOrb(aps,orb,false);
                }
            }
            return;
        }
    
    
    void QMMInterface::addMMAtomtoOrb(ctp::APolarSite * aps,Orbitals &orb, bool with_polarisation){
        vec pos = aps->getPos() *tools::conv::nm2ang;
        double Q = aps->getQ00();
        orb.AddAtom(aps->getName(), pos, Q, true);

        if (_split_dpl) {
            vec tot_dpl=vec(0,0,0);
            if( with_polarisation){
            tot_dpl += aps->getU1();
            }
            if (aps->getRank() > 0) {
                tot_dpl += aps->getQ1();
            }
            // Calculate virtual charge positions
            double a = _dpl_spacing; // this is in nm
            double mag_d = abs(tot_dpl); // this is in e * nm
            if (mag_d > 1e-9) { 
                vec dir_d_0 = tot_dpl.normalize();
                vec dir_d = dir_d_0.normalize();
                vec A = pos + 0.5 * a * dir_d *tools::conv::nm2ang; // converted to AA
                vec B = pos - 0.5 * a * dir_d *tools::conv::nm2ang;
                double qA = mag_d / a;
                double qB = -qA;
                orb.AddAtom("A", A, qA, true);
                orb.AddAtom("B", B, qB, true);
            }
        }
        
        if (aps->getRank() > 1 && _split_dpl) {
            tools::matrix components=aps->getQ2cartesian();
            tools::matrix::eigensystem_t system;
            components.SolveEigensystem(system);
            double a = _dpl_spacing;
            string Atomnameplus[] = {"Xp","Yp","Zp"};
            string Atomnameminus[] = {"Xm","Ym","Zm"};
            for(unsigned i=0;i<3;i++){
            
                double q=system.eigenvalues[i]/(a*a);
                if(q<1e-9){continue;}
                tools::vec vec1=pos+0.5*a*system.eigenvecs[i]*tools::conv::nm2ang;
                tools::vec vec2=pos-0.5*a*system.eigenvecs[i]*tools::conv::nm2ang; 
                orb.AddAtom(Atomnameplus[i], vec1, q, true);
                orb.AddAtom(Atomnameminus[i], vec2, q, true);
            
            }
        }

        return;
    }
    


    }
}
