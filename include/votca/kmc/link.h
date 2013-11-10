/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_KMC_LINK_H
#define	_VOTCA_KMC_LINK_H

#include <votca/tools/vec.h>

namespace votca { namespace kmc {

enum CarrierType{ Electron, Hole, Exciton };    
    
/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separation 
 * 
 */
class Link
{

public:
    /// r2 - r1
    const votca::tools::vec &dr() const { return dr; }    
    
    /// forward and backward rates
    const double &rate12() const { return rate12; }
    const double &rate21() const { return rate21; }
    
private:
    /// forward and backward rates
    double rate12;
    double rate21;
    
    /// r2 - r1
    votca::tools::vec dr;       
};

}} 

#endif // _VOTCA_KMC_LINK_H
