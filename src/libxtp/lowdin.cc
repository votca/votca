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

#include <votca/xtp/lowdin.h>
#include <votca/xtp/aomatrix.h>

namespace votca { namespace xtp {

void Lowdin::EvaluateLowdin(Orbitals& orbitals,const AOBasis &basis,const QMState& state){
    AOOverlap overlap;
    // Fill overlap
    overlap.Fill(basis);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap.Matrix());
    Eigen::MatrixXd Smsqrt=es.operatorSqrt();
    Eigen::MatrixXd prodmat=Smsqrt*orbitals.DensityMatrixFull(state)*Smsqrt;

    int id =0;
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
         PolarSite site=PolarSite(atom);
         site.setCharge(charge);
         orbitals.Multipoles().push_back(site);
         id+=nooffunc;
    }

    return;
}

}}
