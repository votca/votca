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

#ifndef VOTCA_XTP_POPULATIONANALYSIS_H
#define VOTCA_XTP_POPULATIONANALYSIS_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/qmfragment.h>


/**
* \brief Takes a list of atoms, and the corresponding density matrix and puts out a table of wavefunction partial charges
*
* 
* 
*/

namespace votca { namespace xtp {


template <bool T>
class Populationanalysis{
public:
    
    void CalcChargeperAtom(Orbitals& orbitals,const AOBasis &basis,const QMState& state)const{
        Eigen::MatrixXd dmat=orbitals.DensityMatrixFull(state);
        Eigen::VectorXd charges=CalcElecChargeperAtom(dmat,basis);
        if(!state.isTransition()){
            charges+=CalcNucChargeperAtom(orbitals.QMAtoms());
        }

        StaticSegment seg=StaticSegment(orbitals.QMAtoms().getName(),orbitals.QMAtoms().getId());
        for (int i=0;i<orbitals.QMAtoms().size();++i){
             seg.push_back(PolarSite(orbitals.QMAtoms()[i],charges(i)));
        }
        orbitals.Multipoles()=seg;
        return;
    }
    template< typename S >
    Eigen::VectorXd CalcChargeperFragment(const std::vector<QMFragment<S> >& frags, const Eigen::VectorXd& atomcharges){
        Eigen::VectorXd result=Eigen::VectorXd::Zero(frags.size());

        return result;
    }


    private:

        Eigen::VectorXd CalcNucChargeperAtom(const QMMolecule& mol)const{
            Eigen::VectorXd result=Eigen::VectorXd::Zero(mol.size());
            for( int i=0;i<mol.size();i++){
                result(i)=mol[i].getNuccharge();
            }
            return result;
        }

        Eigen::VectorXd CalcElecChargeperAtom(const Eigen::MatrixXd& dmat,const AOBasis &basis)const{
            AOOverlap overlap;
            // Fill overlap
            overlap.Fill(basis);
            Eigen::MatrixXd prodmat;
            if(T){
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap.Matrix());
                Eigen::MatrixXd Smsqrt=es.operatorSqrt();
                prodmat=Smsqrt*dmat*Smsqrt;
            }else{
                prodmat =dmat* overlap.Matrix();
            }
            int noofatoms=basis.getFuncPerAtom().size();
            Eigen::VectorXd charges=Eigen::VectorXd::Zero(noofatoms);
            int start=0;
            for(int i=0;i<charges.size();++i){
                int nofunc=basis.getFuncPerAtom()[i];
                charges(i)=prodmat.diagonal().segment(start,nofunc).sum();
                start+=nofunc;
            }
            return charges;
        }
        

        
        

    
};

typedef Populationanalysis<false> Mulliken;
typedef Populationanalysis<true> Lowdin;



}}

#endif	// VOTCA_XTP_POPULATIONANALYSIS_H
