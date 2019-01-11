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

#ifndef VOTCA_XTP_COUPLINGBASE_H
#define	VOTCA_XTP_COUPLINGBASE_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/logger.h>
#include <boost/format.hpp>
#include <votca/xtp/aomatrix.h>


namespace votca { namespace xtp {


/**
* \brief Base Class to derive DFT and BSE coupling from
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
*/

class CouplingBase
{
public:
    
    
    virtual void CalculateCouplings(const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB, 
                               Orbitals& orbitalsAB)=0;
    
    
    virtual void Initialize(tools::Property&)=0;
    
    virtual void Addoutput(tools::Property & type_summary,const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB)=0;
    
    void setLogger( Logger* pLog ) { _pLog = pLog; }
    
protected:
    Logger *_pLog;
    void CheckAtomCoordinates(const Orbitals& orbitalsA, const Orbitals& orbitalsB, const Orbitals& orbitalsAB);
   
   Eigen::MatrixXd CalculateOverlapMatrix(const Orbitals& orbitalsAB);
  

};

inline Eigen::MatrixXd CouplingBase::CalculateOverlapMatrix(const Orbitals& orbitalsAB){
  BasisSet dftbasisset;
  AOBasis dftbasis;
  dftbasisset.LoadBasisSet(orbitalsAB.getDFTbasisName());
  dftbasis.AOBasisFill(dftbasisset, orbitalsAB.QMAtoms());
  AOOverlap dftAOoverlap;
  dftAOoverlap.Fill(dftbasis);
  Eigen::MatrixXd overlapAB=dftAOoverlap.Matrix();
  return overlapAB;
}

inline void CouplingBase::CheckAtomCoordinates(const Orbitals& orbitalsA, 
                          const Orbitals& orbitalsB, const Orbitals& orbitalsAB){
  const QMMolecule& atomsA=orbitalsA.QMAtoms();
  const QMMolecule& atomsB=orbitalsB.QMAtoms();
  const QMMolecule& atomsAll = orbitalsAB.QMAtoms();
  bool coordinates_agree=true;
  for (int i = 0; i < atomsAll.size(); i++) {
    const QMAtom& dimer = atomsAll[i];
    const QMAtom* monomer = NULL;
    
    if (i < atomsA.size()) {
      monomer = &atomsA[i];
    } else if (i < atomsB.size() + atomsA.size()) {
      monomer = &atomsB[i - atomsA.size()];
    } else {
      // Linker
      XTP_LOG(logERROR, *_pLog) << (boost::format("Neither Monomer A nor Monomer B contains "
              "atom %s on line %u. Hence, this atom is part of a linker.") %dimer.getElement() %(i+1) ).str()<<std::flush;
      continue;
    }

    if (!monomer->getPos().isApprox(dimer.getPos(), 0.001)){
        coordinates_agree=false;
    }
    
    if (monomer->getElement() != dimer.getElement()) {
      throw std::runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
    }
  }
  
  if(!coordinates_agree){
        XTP_LOG(logINFO, *_pLog) << "======WARNING=======\n Coordinates of monomer "
              "and dimer atoms do not agree" << std::flush;
  }
}


}}

#endif	// VOTCA_XTP_COUPLINGBASE_H



