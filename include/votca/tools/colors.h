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

    How to use: std::cout << tools::colors::red << "Hello world";  
*/

namespace Colors 
{
    static const char Reset[]       = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };

    static const char Black[]       = { 0x1b, '[', '0', ';', '3', '0', 'm', 0 };
    static const char Red[]         = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
    static const char Green[]       = { 0x1b, '[', '0', ';', '3', '2', 'm', 0 };
    static const char Yellow[]      = { 0x1b, '[', '0', ';', '3', '3', 'm', 0 };
    static const char Blue[]        = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
    static const char Magenta[]     = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
    static const char Cyan[]        = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
    static const char White[]       = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };

    static const char BoldBlack[]   = { 0x1b, '[', '1', ';', '3', '0', 'm', 0 };
    static const char BoldRed[]     = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
    static const char BoldGreen[]   = { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
    static const char BoldYellow[]  = { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
    static const char BoldBlue[]    = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
    static const char BoldMagenta[] = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
    static const char BoldCyan[]    = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
    static const char BoldWhite[]   = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
        
};

}}

#endif	/* __VOTCA_TOOLS_COLORS_H */

