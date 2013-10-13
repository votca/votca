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

#ifndef _VOTCA_TOOLS_PROPERTY_FORMAT_H
#define	_VOTCA_TOOLS_PROPERTY_FORMAT_H

#include <votca/tools/property.h>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace votca { namespace tools {
      
/**
  * \brief Manipulates the format state of the output stream 
  *
  * Changes the state of the output stream. Property class reads this state  
  * and formats its output according to this state (XML, TXT, T2T, etc)
  */
class PropertyFormat {
    
public:

    /// types of output 
    enum Format{ XML, TXT, LOG, TEX, T2T, HLP };
    
    explicit PropertyFormat(Format fmt, int level = 0, int indent = 0 ) :
                                  _fmt(fmt), 
                                  _level(level),
                                  _indent(indent) {}    

    friend std::ostream& operator << (std::ostream& os, const PropertyFormat& pf)
    {
        os.iword(Property::outputFormat()) = pf._fmt;
        os.iword(Property::outputLevel()) = pf._level;
        os.iword(Property::outputIndent()) = pf._indent;
       return os;
    }

    static void PrintNodeTXT(std::ostream &out, Property &p, const int start_level, int level=0, string prefix="", string offset="");
    static void PrintNodeXML(std::ostream &out, Property &p, const int start_level, int level=0, string prefix="", string offset="");
    static void PrintNodeT2T(std::ostream &out, const string &prefix, Property &p);
    static void PrintNodeLOG(std::ostream &out, const string &prefix, Property &p);
    static void PrintNodeTEX(std::ostream &out, Property &p, const int start_level, int level=0, string prefix="", string offset="");
    static void PrintNodeHLP(std::ostream &out, Property &p, const int start_level, int level=0, string prefix="",  int offset=0);    

private:

    int _fmt; 
    int _level;
    int _indent;
        
};

/**
  * \brief Manipulates the XML-level state of the output stream 
  *
  * forces property object to output nodes starting from a certain level 
  */
class setlevel {
    
public:
    explicit setlevel(int level) : _level(level){}    
    friend std::ostream& operator << (std::ostream& os, const setlevel& pl)
    {
        os.iword(Property::outputLevel()) = pl._level;
        return os;
    }
private:
    int _level;     
};

extern PropertyFormat XML;
extern PropertyFormat TXT;
extern PropertyFormat T2T;
extern PropertyFormat LOG;
extern PropertyFormat TEX;
extern PropertyFormat HLP;


}}


#endif	/* _VOTCA_TOOLS_PROPERTY_FORMAT_H */
