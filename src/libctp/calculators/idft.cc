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

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    ParseOptionsXML( options );
    _orbitalsA.ReadOrbitalsGaussian( _orbitalsA_file.c_str() );
}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {

    string key = "options.idft";
    
    if ( opt->exists(key+".orbitals_A") ) {
        _orbitalsA_file = opt->get(key + ".orbitals_A").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule A orbitals filename is missing.");
    }
    
    /* --- ORBITALS.XML Structure ---
     * <options>
     *   <idft>
     *     <orbitals_A>fort.7</orbitals_A>
     *   </idft>
     * </options>
     */

}

/*
void IDFT::CleanUp() {

}
*/


/*
void IDFT::CalculateJ(QMPair *pair) {

}
*/
    
    
}};