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


#include "idft.h"
#include "votca/ctp/qmcalculator.h"

namespace votca { namespace ctp {
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize( tools::Property* options ) {
        ParseOrbitalsXML( options );
        _orbitalsA.ReadOrbitalsGaussian( "fort.7" );
}

    
void IDFT::ParseOrbitalsXML( tools::Property *opt ) {

    string key = "options.orbitals";
    string orbitalsXML = opt->get(key+".orbitalsXML").as<string> ();
    cout << endl << "... ... Orbital data from " << orbitalsXML << ". ";

    Property alloc;
    load_property_from_xml(alloc, orbitalsXML.c_str());    

    /* --- ORBITALS.XML Structure ---
     */
 
    key = "topology.molecules.molecule";
    list<Property*> mols = alloc.Select(key);
    list<Property*> ::iterator molit;

}

/*
void IDFT::CleanUp() {

}
*/

/*
void IDFT::Initialize(Topology *top, Property *options) {

     cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;

    this->ParseOrbitalsXML(top, options);

    
}
*/

/*
void IDFT::ParseOrbitalsXML(Topology *top, Property *opt) {

}
*/

/*
void IDFT::EvalPair(Topology *top, QMPair *qmpair, int slot) {

}
*/

/*
void IDFT::CalculateJ(QMPair *pair) {

}
*/
    
    
}};