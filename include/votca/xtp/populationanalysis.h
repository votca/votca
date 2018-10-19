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
  
    void Evaluate(Orbitals& orbitals,const AOBasis &basis,const QMState& state){
        AOOverlap overlap;
            // Fill overlap
        overlap.Fill(basis);
        Eigen::MatrixXd prodmat;
        if(T){
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap.Matrix());
            Eigen::MatrixXd Smsqrt=es.operatorSqrt();
            prodmat=Smsqrt*orbitals.DensityMatrixFull(state)*Smsqrt;
        }else{
            prodmat = orbitals.DensityMatrixFull(state)* overlap.Matrix();
        }
        int id =0;
        PolarSegment seg=PolarSegment(orbitals.QMAtoms().getName(),orbitals.QMAtoms().getId());
        for (const QMAtom& atom:orbitals.QMAtoms()){
            double charge=0.0;
             // get element type and determine its nuclear charge
             if (!state.isTransition()){
                charge=atom.getNuccharge();
             }

             int nooffunc=basis.getFuncOfAtom(atom.getAtomID());

             for ( int i = id ; i < id+nooffunc; i++){
                    charge -= prodmat(i,i);
            }
             seg.push_back(PolarSite(atom,charge));
             id+=nooffunc;
        }
        orbitals.Multipoles()=seg;
        return;
    }

    
};

typedef Populationanalysis<false> Mulliken;
typedef Populationanalysis<true> Lowdin;



}}

#endif	// VOTCA_XTP_POPULATIONANALYSIS_H
