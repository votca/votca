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

#ifndef VOTCA_XTP_QMMINTERFACE_H
#define	VOTCA_XTP_QMMINTERFACE_H

#include <votca/tools/elements.h>

namespace votca { 
  
  namespace xtp {

    class QMAtom;
    class APolarSite;
    class PolarSeg;
    class Segment;
    class PolarTop;
    class Orbitals;

// ========================================================================== //
// QM-MM INTERFACE CLASS - CONVERTS BETWEEN QMATOMS <> POLAR OBJECTS          //
// ========================================================================== //
    
class QMInterface
{
public:
    
    // CONVERSION QM -> MM
    xtp::APolarSite *Convert(QMAtom *atm, int id = -1);
    
    xtp::PolarSeg Convert(std::vector<QMAtom*> &atms);
    
    std::vector<QMAtom *> Convert( std::vector<xtp::Segment* > segments);
    
    void GenerateQMAtomsFromPolarSegs(xtp::PolarTop *ptop, Orbitals &orb);
    std::vector<std::shared_ptr<xtp::PolarSeg> > GenerateMultipoleList(xtp::PolarTop *ptop  );
    void Orbitals2Segment(xtp::Segment& segment, const Orbitals& orbitals);
     
private:
    void addMMAtomtoOrb(xtp::APolarSite * aps,Orbitals &orb, bool with_polarisation);
    // Allocates polarizabilities in A**3 to element types
    tools::Elements _element;
   
};




}}

#endif // VOTCA_XTP_QMMINTERFACE_H
