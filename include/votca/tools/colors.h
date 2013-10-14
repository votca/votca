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


#ifndef __VOTCA_TOOLS_COLORS_H
#define	__VOTCA_TOOLS_COLORS_H

namespace votca { namespace tools {

/**
    \brief class to store color codes     

    Use: std::cout << tools::colors::red;  
*/

struct colors 
{
        /// colors USE: std::cout << tools::colors::red
        struct colors {
             char normal[]      = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
             char black[]       = { 0x1b, '[', '0', ';', '3', '0', 'm', 0 };
             char black_bold[]  = { 0x1b, '[', '1', ';', '3', '0', 'm', 0 };
             char blue[]        = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
             char blue_bold[]   = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
             char blue_dark[]   = { 0x1b, '[', '2', ';', '3', '4', 'm', 0 };

             char red[]         = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
             char green[]       = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
             char red[]         = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
             char red[]         = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
};

}}

#endif	/* __VOTCA_TOOLS_COLORS_H */

