/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
// 
// File:   boltzmanninversion.h
// Author: victor
//
// Created on 4. Juni 2008, 17:39
//

#ifndef _BONDEDSTATISTICS_H
#define	_BONDEDSTATISTICS_H

#include "cgobserver.h"
#include <votca/tools/datacollection.h>

class BondedStatistics
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom = 0);
    void EndCG();
    
    void EvalConfiguration(Topology *conf, Topology *conf_atom = 0);
    
    DataCollection<double> &BondedValues() { return _bonded_values; }

protected:
    DataCollection<double> _bonded_values;
};

#endif	/* _BOLZMANNINVERSION_H */

