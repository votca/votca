/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_BONDEDSTATISTICS_H
#define	VOTCA_CSG_BONDEDSTATISTICS_H

#include <votca/csg/cgobserver.h>
#include <votca/tools/datacollection.h>

namespace TOOLS = votca::tools;

namespace votca {
  namespace csg {
/**
 * \brief Class calculates data associated with bond interactions
 *
 * This class will essentially calculate and store information related to each
 * type of interaction. For instance for the IBond which is the interaction
 * between two beads it will calculate and store the distance between the two
 * beads involved in the interaction. It will calculate a similar metric for all
 * other interactions such as IAngle, IDihedral etc...
 **/
class BondedStatistics
    : public votca::csg::CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom = 0);
    void EndCG();
    
    void EvalConfiguration(Topology *conf, Topology *conf_atom = 0);
    
    TOOLS::DataCollection<double> &BondedValues() { return _bonded_values; }

protected:
    TOOLS::DataCollection<double> _bonded_values;
};
}}
#endif	// VOTCA_CSG_BOLZMANNINVERSION_H 

