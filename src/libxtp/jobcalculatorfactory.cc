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


#include <votca/xtp/jobcalculatorfactory.h>
#include "votca_config.h"


#include "jobcalculators/idft.h"
#include "jobcalculators/edft.h"
#include "jobcalculators/xqmultipole.h"
#include "jobcalculators/qmmm.h"
#include "jobcalculators/ewald.h"
#include "jobcalculators/egwbse.h"
#include "jobcalculators/igwbse.h"
#include "jobcalculators/qmape.h"
#include "jobcalculators/iexcitoncl.h"


namespace votca { namespace xtp {

void JobCalculatorfactory::RegisterAll(void)
{	
        JobCalculators().Register<EDFT>                ("edft");
        JobCalculators().Register<IDFT>                ("idft");
        JobCalculators().Register<XQMP>                ("xqmultipole");
        JobCalculators().Register<QMMM>                ("qmmm");
        JobCalculators().Register< Ewald<Ewald3D2D> >  ("ewald2d");
        JobCalculators().Register< Ewald<Ewald3D3D> >  ("ewald3d");
        JobCalculators().Register< Ewald<PEwald3D3D> > ("pewald3d");
        JobCalculators().Register<EGWBSE>              ("egwbse");
        JobCalculators().Register<IGWBSE>              ("igwbse");
        JobCalculators().Register<QMAPE>               ("qmape");
        JobCalculators().Register<IEXCITON>            ("iexcitoncl");
}

}}
