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
#include "calculators/panalyze.h"
#include "calculators/eimport.h"
#include "calculators/pairdump.h"
#include "calculators/iimport.h"
#include "calculators/ianalyze.h"
#include "calculators/profile.h"
#include "calculators/xmultipole2.h"
#include "calculators/velocity.h"
#include "calculators/zmultipole.h"
#include "calculators/qmultipole.h"
#include "calculators/jobwriter.h"
#include "calculators/ewdbgpol.h"
#include "calculators/cgpolar.h"
#include "calculators/vaverage.h"




namespace votca { namespace xtp {

void Calculatorfactory::RegisterAll(void)
{	
        Calculators().Register<Sandbox>             ("sandbox");
        Calculators().Register<Neighborlist>        ("neighborlist");
        Calculators().Register<StateServer>         ("stateserver");
        Calculators().Register<TDump>               ("tdump");
        Calculators().Register<IZindo>              ("izindo");
        Calculators().Register<Rates>               ("rates");
        Calculators().Register<EInternal>           ("einternal");
        Calculators().Register<EOutersphere>        ("eoutersphere");
        Calculators().Register<EMultipole>          ("emultipole");
        Calculators().Register<EAnalyze>            ("eanalyze");
        Calculators().Register<PAnalyze>            ("panalyze");
        Calculators().Register<EImport>             ("eimport");
        Calculators().Register<PairDump>            ("pairdump");        
        Calculators().Register<IImport>             ("iimport");
        Calculators().Register<IAnalyze>            ("ianalyze");
        Calculators().Register<Profile>             ("profile");
        Calculators().Register<XMP>                 ("xmultipole");
        Calculators().Register<Velocity>            ("velocity");
        Calculators().Register<ZMultipole>          ("zmultipole");
        Calculators().Register<QMultipole>          ("qmultipole");
        Calculators().Register<JobWriter>           ("jobwriter");
        Calculators().Register<EwaldBgPolarizer>    ("ewdbgpol");
        Calculators().Register<CgPolar>             ("cgpolar");
        Calculators().Register<VAverage>            ("vaverage");
}

}}
