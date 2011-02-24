/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "trajectoryreader.h"
#include "modules/io/esptrajectoryreader.h"
#include "modules/io/lammpsreader.h"

#ifdef GMX
#include "modules/io/gmxtrajectoryreader.h"
#endif

namespace votca { namespace csg {

void TrajectoryReader::RegisterPlugins(void)
{
    TrjReaderFactory().Register<ESPTrajectoryReader>("esp");
    TrjReaderFactory().Register<LAMMPSReader>("dump");
#ifdef GMX
    TrjReaderFactory().Register<GMXTrajectoryReader>("trr");
    TrjReaderFactory().Register<GMXTrajectoryReader>("xtc");
    TrjReaderFactory().Register<GMXTrajectoryReader>("gro");
    TrjReaderFactory().Register<GMXTrajectoryReader>("pdb");
#endif
}

}}
