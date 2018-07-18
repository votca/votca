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

#ifndef __QMINTERFACE__H
#define	__QMINTERFACE__H




#include <votca/ctp/apolarsite.h>
#include <votca/xtp/qmatom.h>
#include <votca/ctp/polarseg.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/polartop.h>
// add gwbse header for excited state support
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/gdma.h>


namespace votca { namespace xtp {

    
// ========================================================================== //
// QM-MM INTERFACE CLASS - CONVERTS BETWEEN QMATOMS <> POLAR OBJECTS          //
// ========================================================================== //
    
class QMInterface
{
public:
    
    QMInterface() { _polar_table = ctp::POLAR_TABLE(); };
   ~QMInterface() {};
    
    // CONVERSION QM -> MM
    ctp::APolarSite *Convert(QMAtom *atm, int id = -1);
    
    ctp::PolarSeg Convert(std::vector<QMAtom*> &atms);
    
    void setMultipoleSplitting(bool split_dpl, double dpl_spacing){
        _split_dpl=split_dpl;
        _dpl_spacing=dpl_spacing;
    }
    
    std::vector<QMAtom *> Convert( std::vector<ctp::Segment* > segments);
    
    void GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop, Orbitals &orb);
    std::vector<ctp::PolarSeg*> GenerateMultipoleList(ctp::PolarTop *ptop  );
    void Orbitals2Segment(ctp::Segment* _segment, const Orbitals& _orbitals);
    
     
private:
    void addMMAtomtoOrb(ctp::APolarSite * aps,Orbitals &orb, bool with_polarisation);
    // Allocates polarizabilities in A**3 to element types
    std::map<std::string,double> _polar_table;
    bool _split_dpl;
    double _dpl_spacing;
};




}}

#endif
