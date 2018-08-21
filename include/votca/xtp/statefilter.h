/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_STATEFILTER_H
#define _VOTCA_XTP_STATEFILTER_H

#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>


namespace votca {
namespace xtp {
/**
 *  \brief  Filters from a spectrum of states the state, which fullfills certain criteria
 *
 *
 */

class Statefilter {

public:
    void Initialize(tools::Property& options);
    
    void setType(const std::string& type){_type=type;}
    
    void Filter(const Orbitals& orbital);
    
    int getStateIndex(){return _state_index;}// zero indexed;
    
private:
 
std::string type;
int _initial_state_index;
int _state_index;    
ctp::Logger *_log;
 




};
}
}

#endif /* _VOTCA_XTP_STATEFILTER_H */
