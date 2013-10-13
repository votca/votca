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

#include <votca/tools/propertyformat.h>

namespace votca { namespace tools {

PropertyFormat XML(PropertyFormat::XML);
PropertyFormat TXT(PropertyFormat::TXT);
PropertyFormat T2T(PropertyFormat::T2T);
PropertyFormat LOG(PropertyFormat::LOG);
PropertyFormat TEX(PropertyFormat::TEX);
PropertyFormat HLP(PropertyFormat::HLP);

void PropertyFormat::PrintNodeTXT(std::ostream &out, Property &p, const int start_level, int level, string prefix, string offset)
{
    
    list<Property>::iterator iter;
        
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

void PropertyFormat::PrintNodeXML(std::ostream &out, Property &p, const int start_level, int level, string prefix,  string offset)
{
    
    list<Property>::iterator iter;       
    Property::AttributeIterator ia;
    bool _endl = true;
    
        // print starting only from the start_level (the first node (level 0) can be <> </>)
        if ( level >= start_level )  {
            // print the node name
            out << offset << "<" << p.name() ;
            // print the node attributes 
            for(ia = p.firstAttribute(); ia!=p.lastAttribute(); ++ia) 
                out << " " << ia->first << "=\"" << ia->second << "\"" ;
            out << ">";
            
            // print node value if it is not empty
            bool has_value = ( (p.value()).find_first_not_of("\t\n ") != std::string::npos );
            if( has_value ) { out << p.value(); _endl = false; }
            
            // check if we need the end of the line or not
            if( !has_value &&  p.HasChilds() ) out << endl;
            if( !has_value &&  !p.HasChilds() ) _endl = false;
        }
    
        // continue iteratively through the rest of the nodes
        for(iter = p.begin(); iter!=p.end(); ++iter) {
            level++; 
            if (  level > start_level ) offset += "\t"; 
            PrintNodeXML(out, (*iter), start_level, level, (*iter).name(), offset);
            if (  level > start_level ) offset.resize(offset.size()-1); 
            level--;           
        }
        
        if ( level >= start_level ) {
            if ( _endl ) {
                out << offset << "</" << p.name() << ">"  << endl;
            } else {
                out << "</" << p.name() << ">"  << endl;
            }
        } 
}

void PropertyFormat::PrintNodeT2T(std::ostream &out, const string &prefix, Property &p) {
    out << "T2T format is not implemented\n";
}

void PropertyFormat::PrintNodeLOG(std::ostream &out, const string &prefix, Property &p) {
    out << "LOG format is not implemented\n";
}
    
void PropertyFormat::PrintNodeTEX(std::ostream &out, Property &p, const int start_level, int level, string prefix,  string offset) {

    
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
            head_name = p.name();
            if ( p.hasAttribute("label") ) _label = p.getAttribute<string>("label");
            if ( p.hasAttribute("section") ) _section = p.getAttribute<string>("section");
            if ( p.hasAttribute("help") ) _help = p.getAttribute<string>("help");
            out << boost::format(header_format) % head_name % _label % _help;     
            prefix = p.name();
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
            PrintNodeTEX(out, (*iter), start_level, level, (*iter).name(), offset);
            level--;
        } else {
            level++;
            PrintNodeTEX(out, (*iter), start_level, level, prefix + "." + (*iter).name(), offset);
            level--;
        }
    }        

    // if this is the head node, print the footer
    if ( level == start_level )  out << boost::format(footer_format) % _section % head_name;
}

void PropertyFormat::PrintNodeHLP(std::ostream &out, Property &p, const int start_level, int level, string prefix,  int offset) {

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
            PropertyFormat::PrintNodeTXT(out, p, _level);
            case PropertyFormat::XML:
            PropertyFormat::PrintNodeXML(out, p, _level);
            break;
        case PropertyFormat::TXT:
            PropertyFormat::PrintNodeTXT(out, p, _level, 0, "", sindent);
            break;
        case PropertyFormat::T2T:
            PropertyFormat::PrintNodeT2T(out, "", p);
            break;
        case PropertyFormat::LOG:
            PropertyFormat::PrintNodeLOG(out, "", p);
            break;
        case PropertyFormat::TEX:            
            PropertyFormat::PrintNodeTEX(out, p, _level);
            break;
        case PropertyFormat::HLP:            
            PropertyFormat::PrintNodeHLP(out, p, _level);
            break;
        }
        //out << endl;
    }

    return out;
};


}}
