/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#ifndef _CALC_DFT_ENERGIES_H
#define	_CALC_DFT_ENERGIES_H

#include <votca/ctp/segment.h>
#include <votca/ctp/parallelsitecalc.h>

namespace votca { namespace ctp {

/**
* \brief Site energies and orbitals for QM pairs
*
* QM orbitals and energies for all molecules
* Requires first-principles package, i.e. GAUSSIAN installation
*
* Callname: idft
*/

class EDFT : public ParallelSiteCalculator
{
public:

    EDFT() {};
   ~EDFT() {};

    string  Identify() { return "EDFT"; }
    void    Initialize(Topology *top, Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    void    EvalSite(Topology *top, Segment *seg, int slot);

    void    CleanUp();


private:

    bool        _maverick;

};

void EDFT::CleanUp() {

}

void EDFT::Initialize(Topology *top, Property *options) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;
    cout << endl <<  "... ... Reading the input ";
    
    /* ---- OPTIONS.XML Structure -----
     *
     * <edft>
     *
     *      <package>gaussian.xml</package>
     *
     * </edft>
     *
     */

    //this->ParseOrbitalsXML(top, options);

}


void EDFT::ParseOrbitalsXML(Topology *top, Property *opt) {

    string key = "options.edft";
    string program = opt->get(key+".method").as<string> ();
    cout << endl << "... ... Using " << program << ". ";

    Property alloc;
    load_property_from_xml(alloc, program.c_str());    

}


void EDFT::EvalSite(Topology *top, Segment *seg, int slot) {

    exit(0);
    this->LockCout();
    //cout << "\r... ... Evaluating site " << seg->getId()+1 << flush;   
    //string SiteName = seg->getName();
    //cout << "... ... ... Site name " << SiteName << endl; 
    //cout << "... ... ... slot " << slot << endl;    
    this->UnlockCout();

}

}}

#endif	/* _CALC_DFT_ENERGIES_H */
