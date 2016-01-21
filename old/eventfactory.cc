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


#include <votca/kmc/eventfactory.h>
#include "votca_config.h"
#include "events/electron_transfer.h"
#include "events/hole_transfer.h"
#include "events/electron_injection.h"
#include "events/hole_injection.h"
#include "events/electron_collection.h"
#include "events/hole_collection.h"
#include "events/recombination.h"

namespace votca { namespace kmc {

void EventFactory::RegisterAll(void)
{
<<<<<<< local
        //Events().Register< ElectronTransfer >         ( _ElectronTransfer );
        //Events().Register< HoleTransfer >         ( _HoleTransfer );
        //Events().Register< ElectronInjection >         ( _ElectronInjection );
        //Events().Register< HoleInjection >         ( _HoleInjection );
        //Events().Register< ElectronCollection >         ( _ElectronCollection );
        //Events().Register< HoleCollection >         ( _HoleCollection );
        //Events().Register< Recombination >         ( _Recombination );
=======
        Events().Register< ElectronTransfer >         ( _ElectronTransfer );
        Events().Register< HoleTransfer >         ( _HoleTransfer );
        Events().Register< ElectronInjection >         ( _ElectronInjection );
        Events().Register< HoleInjection >         ( _HoleInjection );
        Events().Register< ElectronCollection >         ( _ElectronCollection );
        Events().Register< HoleCollection >         ( _HoleCollection );
        Events().Register< Recombination >         ( _Recombination );
>>>>>>> other
}

}}

