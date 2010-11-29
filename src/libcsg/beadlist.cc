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

#include "beadlist.h"
#include "topology.h"
#include <votca/tools/tokenizer.h>

namespace votca { namespace csg {

int BeadList::Generate(Topology &top, const string &select)
{
    BeadContainer::iterator iter;
    _topology = &top;
    
    for(iter=top.Beads().begin(); iter!=top.Beads().end();++iter) {
        if(wildcmp(select.c_str(), (*iter)->getType()->getName().c_str())) {
            push_back(*iter);
        }
    }
    return size();
}

}}