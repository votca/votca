/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#include <votca/xtp/calculatorfactory.h>

//#include "calculators/sandbox.h"
#include "calculators/neighborlist.h"
/*
 #include "calculators/stateserver.h"
#include "calculators/rates.h"
#include "calculators/einternal.h"
#include "calculators/eanalyze.h"
#include "calculators/panalyze.h"
#include "calculators/eimport.h"
#include "calculators/iimport.h"
#include "calculators/ianalyze.h"
*/




namespace votca { namespace xtp {
    namespace CTP = votca::ctp;
    
void XCalculatorfactory::RegisterAll(void)
{	
        //XCalculators().Register<votca::ctp::Sandbox>             ("sandbox");
	XCalculators().Register<XNeighborlist>       ("xneighborlist");
        /* XCalculators().Register<XStateServer>        ("xstateserver");
        XCalculators().Register<Rates>               ("rates");
        XCalculators().Register<XEInternal>          ("xeinternal");
        XCalculators().Register<XEAnalyze>           ("xeanalyze");
        XCalculators().Register<PAnalyze>            ("panalyze");
        XCalculators().Register<XEImport>            ("xeimport");
        XCalculators().Register<XIImport>            ("xiimport");
        XCalculators().Register<XIAnalyze>           ("ianalyze");
        */

}

}}
