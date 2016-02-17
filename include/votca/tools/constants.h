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

#ifndef __VOTCA_TOOLS_CONSTANTS_H
#define	__VOTCA_TOOLS_CONSTANTS_H

#include <boost/math/constants/constants.hpp>

namespace votca { namespace tools {
    
    
    namespace conv{
        
        // mathematical constants 
        
    const double Pi = boost::math::constants::pi<double>();
        // natural constants
        
    const double kB = 8.617332478E-5; // eV/K
    const double hbar = 6.5821192815E-16; // eV*s
        
        //length conversions
   //votca xtp-uses for any conversions the following scheme unitA2unitB 
    const double bohr2nm =0.052917721092;
    const double nm2bohr =18.897259886;
    const double ang2bohr = 1.8897259886; 
    const double bohr2ang =1.0/1.8897259886; 
    const double nm2ang=10.0;
    const double ang2nm=0.1;
    
    
    const double ryd2ev= 13.60569253;
    const double ha2ev=  27.21138602;
    
    //ewald internal to eV conversion
    const double int2eV = 1/(4*Pi*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    
    // energy conversions
    static const double rydtoev = 13.6058;
    
    // floats for gwbse
    
    const float ryd2ev_f= 13.60569253;
    
    }

}}

#endif	/* CONVERSIONFACTORS */

