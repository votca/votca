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

class PropertyFormat;   
extern PropertyFormat XML;
extern PropertyFormat TXT;
extern PropertyFormat T2T;
extern PropertyFormat LOG;
extern PropertyFormat TEX;
extern PropertyFormat HLP;
    
    
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

void PrintNodeTXT(std::ostream &out, Property &p, const int start_level, int level, string prefix, string offset)
{
    
    list<Property>::iterator iter;
    
    /* int indent = 10;
    string fmt = "%|" + boost::lexical_cast<string>(indent) +"t|%1%%|15t|%2%\n";
    out << boost::format(fmt) % "A" % "B"; */
    
    if((p.value() != "") || p.HasChilds() ) {
        
        if ( level >= start_level ) {
                if((p.value()).find_first_not_of("\t\n ") != std::string::npos)   
                        out << offset << prefix << " = " << p.value() << endl;
        } else {
            prefix="";
        }
    }
        
    for(iter = p.begin(); iter!=p.end(); ++iter) {
        if(prefix=="") {
            level++;
            PrintNodeTXT(out, (*iter), start_level, level, prefix + (*iter).name(), offset );
            level--;
        } else {
            level++;
            PrintNodeTXT(out, (*iter), start_level, level, prefix + "." + (*iter).name(), offset );
            level--;
        }
    }
}

void PrintNodeXML(std::ostream &out, Property &p, const int start_level, int level, string prefix,  string offset)
{
    
    list<Property>::iterator iter;       
    map<string,string>::iterator ia;
    bool _endl = true;
    
        // print starting only from the start_level (the first node (level 0) can be <> </>)
        if ( level >= start_level )  {
            // print the node name
            out << offset << "<" << p._name ;
            // print the node attributes 
            for(ia = p._attributes.begin(); ia!=p._attributes.end(); ++ia) 
                out << " " << ia->first << "=\"" << ia->second << "\"" ;
            out << ">";
            
            // print node value if it is not empty
            bool has_value = ( (p._value).find_first_not_of("\t\n ") != std::string::npos );
            if( has_value ) { out << p._value; _endl = false; }
            
            // check if we need the end of the line or not
            if( !has_value &&  p.HasChilds() ) out << endl;
            if( !has_value &&  !p.HasChilds() ) _endl = false;
        }
    
        // continue iteratively through the rest of the nodes
        for(iter = p._properties.begin(); iter!=p._properties.end(); ++iter) {
            level++; 
            if (  level > start_level ) offset += "\t"; 
            PrintNodeXML(out, (*iter), start_level, level, (*iter)._name, offset);
            if (  level > start_level ) offset.resize(offset.size()-1); 
            level--;           
        }
        
        if ( level >= start_level ) {
            if ( _endl ) {
                out << offset << "</" << p._name << ">"  << endl;
            } else {
                out << "</" << p._name << ">"  << endl;
            }
        } 
}

void PrintNodeT2T(std::ostream &out, const string &prefix, Property &p) {
    out << "T2T format is not implemented";
}

void PrintNodeLOG(std::ostream &out, const string &prefix, Property &p) {
    out << "LOG format is not implemented";
}
    
void PrintNodeTEX(std::ostream &out, Property &p, const int start_level, int level, string prefix,  string offset) {

    list<Property>::iterator iter;       
    string head_name;
    string _label(""); // reference of the xml file in the manual
    string _section(""); // reference of the description section in the manual
    string _help("");
    string _default(""); // default value if supplied
    string _unit(""); //unit, if supplied
 

    string header_format("\\subsection{%1%}\n"
                         "\\label{%2%}\n%3%\n"
                         "\\rowcolors{1}{invisiblegray}{white}\n"
                         "{\\small\n "
                         "\\begin{longtable}{m{3cm}|m{2cm}|m{1cm}|m{8cm}}\n"
                         " option & default & unit & description\\\\\n\\hline\n");    
    
    string footer_format("\\end{longtable}\n}\n"
                         "\\noindent Return to the description of \\slink{%1%}{\\texttt{%2%}}.\n");
    
    string body_format(" \\hspace{%1%pt}\\hypertarget{%2%}{%3%} & %4% & %5% & %6% \\\\\n");
    
    // if this is the head node, print the header
    if ( level == start_level )  {
            head_name = p._name;
            if ( p.hasAttribute("label") ) _label = p.getAttribute<string>("label");
            if ( p.hasAttribute("section") ) _section = p.getAttribute<string>("section");
            if ( p.hasAttribute("help") ) _help = p.getAttribute<string>("help");
            out << boost::format(header_format) % head_name % _label % _help;     
            prefix = p._name;
     } 
    
    
    if ( level > start_level ) {
    
        // if this node has children or a value or is not the first, start recursive printing
        if( ( p.value() != "" || p.HasChilds() ) && level > -1) {
            string _tex_name = boost::replace_all_copy( p.name(), "_", "\\_" );
            
            if ( p.hasAttribute("default") ) _default = p.getAttribute<string>("default");
            if ( p.hasAttribute("unit") ) _unit = p.getAttribute<string>("unit");
            if ( p.hasAttribute("help") ) _help = p.getAttribute<string>("help");
            
            out << boost::format(body_format) % int((level-start_level-1)*10) 
                    % prefix 
                    % _tex_name 
                    % _default
                    % _unit
                    % _help;
        }
    } 
   
    // continue iteratively through the rest of the nodes
    for(iter = p.begin(); iter != p.end(); ++iter) {
        if(prefix=="") {
            level++;
            PrintNodeTEX(out, (*iter), start_level, level, (*iter)._name, offset);
            level--;
        } else {
            level++;
            PrintNodeTEX(out, (*iter), start_level, level, prefix + "." + (*iter)._name, offset);
            level--;
        }
    }        

    // if this is the head node, print the footer
    if ( level == start_level )  out << boost::format(footer_format) % _section % head_name;

}

void PrintNodeHLP(std::ostream &out, Property &p, const int start_level, int level=0, string prefix="",  int offset=0) {

    list<Property>::iterator iter;       
    string head_name;
    string _help("");
    string _unit("");
    string _default("");
    string fmt("t|%1%%|15t|%2%%|40t|%3%%|55t|%4%\n");
    
    // if this is the head node, print the header
    if ( level == start_level ) {
            head_name = p.name();
            if ( p.hasAttribute("help") ) {
                if  ( p.hasAttribute("help") ) _help = p.getAttribute<string>("help");           
                out << boost::format(" %1%: %|18t| %2%\n") % head_name % _help;
            }
            offset=0;
            //out << boost::format(fmt) % "option" % "def" % "[un]" % "description";
    } 
    
    if ( level > start_level ) {
            
            string ofmt;
            ofmt = "%|" + boost::lexical_cast<string>(offset) + fmt;
            //cout << ofmt << " " << fmt << endl;
            
            if  ( p.hasAttribute("unit") ) _unit = p.getAttribute<string>("unit");
            if  ( p.hasAttribute("default") ) _default = p.getAttribute<string>("default") ;
            if  ( p.hasAttribute("help") ) _help = p.getAttribute<string>("help") ;
            if ( !_unit.empty() ) _unit = "[" + _unit + "]";
            if ( !_default.empty() ) _default = "(" + _default + ")";
            
            out << boost::format(ofmt)
                    % p.name() 
                    % _default
                    % _unit
                    % _help;
    }
    
    for(iter = p.begin(); iter != p.end(); ++iter) {
        if(prefix=="") {
            offset += 2; level++;
            PrintNodeHLP(out, (*iter), start_level, level, prefix + (*iter).name(), offset);
            offset -= 2; level--;
        } else {
            offset += 2; level++;
            PrintNodeHLP(out, (*iter), start_level, level, prefix + "." + (*iter).name(), offset);
            offset -= 2; level--;
        }
    }        

    // if this is the head node, print the footer
    if ( level == start_level ) {
     }
}

std::ostream &operator<<(std::ostream &out, Property& p)
{
    if (!out.good())
        return out;

    std::ostream::sentry sentry(out);

    if(sentry)
    {
        // level from which to start node output
        int _level = out.iword(p.outputLevel());
        int _indent = out.iword(p.outputIndent());
        string sindent;
        for (int i = 0; i < _indent; i++ ) sindent +=  " ";
            
        switch(out.iword(p.outputFormat()))
        {
        default:
            PrintNodeTXT(out, p, _level);
            case PropertyFormat::XML:
            PrintNodeXML(out, p, _level);
            break;
        case PropertyFormat::TXT:
            PrintNodeTXT(out, p, _level, 0, "", sindent);
            break;
        case PropertyFormat::T2T:
            PrintNodeT2T(out, "", p);
            break;
        case PropertyFormat::LOG:
            PrintNodeLOG(out, "", p);
            break;
        case PropertyFormat::TEX:            
            PrintNodeTEX(out, p, _level);
            break;
        case PropertyFormat::HLP:            
            PrintNodeHLP(out, p, _level);
            break;
        }
        //out << endl;
    }

    return out;
};

}}


#endif	/* _VOTCA_TOOLS_PROPERTY_FORMAT_H */
