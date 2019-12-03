/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "modules/io/groreader.h"
#include "modules/io/lammpsdatareader.h"
#include "modules/io/lammpsdumpreader.h"
#include "modules/io/xmltopologyreader.h"
#include <votca/csg/topologyreader.h>
#include <votca/csg/xyzreader.h>
#ifdef GMX_DOUBLE
#include "modules/io/gmxtopologyreader.h"
#endif
#include "modules/io/dlpolytopologyreader.h"
#include "modules/io/pdbreader.h"

namespace votca {
namespace csg {

void TopologyReader::RegisterPlugins(void) {
  TopReaderFactory().Register<XMLTopologyReader>("xml");
  TopReaderFactory().Register<LAMMPSDumpReader>("dump");
  TopReaderFactory().Register<LAMMPSDataReader>("data");
  TopReaderFactory().Register<XYZReader>("xyz");
  TopReaderFactory().Register<GROReader>("gro");
#ifdef GMX_DOUBLE
  TopReaderFactory().Register<GMXTopologyReader>("tpr");
#endif
  TopReaderFactory().Register<PDBReader>("pdb");
  TopReaderFactory().Register<DLPOLYTopologyReader>("dlpf");
}

}  // namespace csg
}  // namespace votca
