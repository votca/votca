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


#include <votca/ctp/calculatorfactory.h>
#include "votca_config.h"

#include "calculators/sandbox.h"
#include "calculators/neighborlist.h"
#include "calculators/stateserver.h"
#include "calculators/tdump.h"
#include "calculators/rates.h"
#include "calculators/izindo.h"
#include "calculators/einternal.h"
#include "calculators/eoutersphere.h"
#include "calculators/emultipole.h"
#include "calculators/emultipole_stdal.h"
#include "calculators/eanalyze.h"




namespace votca { namespace ctp {

void Calculatorfactory::RegisterAll(void)
{	
        Calculators().Register<Sandbox>("sandbox"); // Test calculator
        Calculators().Register<Neighborlist>("neighborlist");
        Calculators().Register<StateServer>("stateserver");
        Calculators().Register<TDump>("tdump");
        Calculators().Register<IZindo>("izindo");
        Calculators().Register<Rates>("rates");
        Calculators().Register<EInternal>("einternal");
        Calculators().Register<EOutersphere>("eoutersphere");
        Calculators().Register<EMultipole>("emultipole");
        Calculators().Register<EMultipole_StdAl>("emultipole2");
        Calculators().Register<EAnalyze>("eanalyze");
}

}}
