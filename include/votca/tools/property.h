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

#ifndef _VOTCA_TOOLS_PROPERTY_H
#define	_VOTCA_TOOLS_PROPERTY_H

#include <map>
#include <string>
#include <iostream>
#include <list>
#include <stdexcept>
#include "lexical_cast.h"
#include <boost/algorithm/string/trim.hpp>
#include <stdlib.h>

#include "vec.h"
#include "tokenizer.h"

namespace votca { namespace tools {
    
/**
 * \brief class to manage program options with xml serialization functionality
 *
 * This class stores tags and content in a hierarchical (tree) structure similar
 * to the one used in the XML format. The structure can be either filled manually
 * or read in from an XML file using load_property_from_xml.
 * The supported XML constructs are TAGS, ATTRIBUTES, and CONTENT:
 * <tag attribute_name="attribute_value">
 *      content
 * </tag> 
 * The property object can be output to an ostream using format modifiers:
 * cout << XML << property;
 * Supported formats are XML, TXT, TEX, HLP
 */
class Property {
    
    /// \brief outputs the property to the ostream
    friend std::ostream &operator<<(std::ostream &out, Property& p);
    /// \brief output the content in the xml format
    friend void PrintNodeXML(std::ostream &out, Property &p, const int start_level=1, int level=0, string prefix="", string offset="");
    /// \brief output the content in the text format
    friend void PrintNodeTXT(std::ostream &out, Property &p, const int start_level=1, int level=0, string prefix="", string offset="");
    /// \brief output the content in the tex format
    friend void PrintNodeTEX(std::ostream &out, Property &p, const int start_level=1, int level=0, string prefix="", string offset="");

   
public:
    Property() : _path("") {}
    
    Property(const string &name, const string &value, const string &path) 
        : _name(name), _value(value), _path(path) {}
    
    /**
     * \brief add a new property to structure
     * @param key identifier
     * @param value value
     * @return reference to the created Property object
     */
    Property &add(const string &key, const string &value);
    /**
     * \brief set value of existing property
     * @param key identifier
     * @param value value
     * @return reference to the created Property object
     */
    Property &set(const string &key, const string &value);
    
    /**
     * \brief get existing property
     * @param key identifier
     * @return Reference to property object
     *
     * This function tries to find a property specified by key separated
     * by "." to step down hierarchy. If the property is not
     * found a runtime_exception is thrown.
     */
    Property &get(const string &key);

    /**
     * \brief check weather property exists
     * @param key identifier
     * @return true or false
     */
    bool exists(const string &key);
    
    /**
     * \brief select property based on a filter
     * @param filter
     * @return list of pointers to property objects
     *
     * returns a list of properties that match the key criteria including
     * wildcards "*" and "?". Example: "base.item*.value"
    */
    std::list<Property *> Select(const string &filter);
    
    /**
     * \brief reference to value of property
     * @return string content
     */
    string &value() { return _value; }
    /**
     * \brief name of property
     * @return name
     */
    string name() { return _name; }
    /**
     * \brief full path of property (including parents)
     * @return path
     *
     * e.g. cg.inverse.value
     */
    string path() { return _path; }
    /**
     * \brief return value as type
     *
     * returns the value after type conversion, e.g.
     * p.as<int>() returns an integer
     */
    template<typename T>
    T as() const;
    /**
     * \brief does the property has childs?
     * \return true or false
     */
    bool HasChilds() { return !_map.empty(); }
    
    /// iterator to iterate over properties
    typedef list<Property>::iterator iterator;  
    /// \brief iterator to first child property
    iterator begin() { return _properties.begin(); }
    /// \brief end iterator for child properties
    iterator end() { return _properties.end(); }
    /// \brief number of child properties
    list<Property>::size_type size() { return _properties.size(); }

    // throw error and comment (with filename+code line)
    void throwRuntimeError(string message);
    /**
     * \brief return attribute as type
     *
     * returns an attribute after type conversion, e.g.
     * p.getAttribute<int>() returns an integer
     */
    template<typename T>
    T getAttribute(const string &attribute);
    /**
     * \brief return attribute as type
     *
     * returns an attribute after type conversion, e.g.
     * p.getAttribute<int>() returns an integer
     */
    template<typename T>
    T getAttribute( std::map<string,string>::iterator it);
    /**
     * \brief set an attribute
     */
    template<typename T>
    void setAttribute(const string &attribute, const T &value);
    /**
     * \brief return true if a node has attributes
     */
    bool hasAttributes() { return _attributes.size() > 0; }
    /**
     * \brief return true if an attribute exists
     */
    bool hasAttribute(const string &attribute);
    /**
     * \brief returns and iterator to an attribute
     */    
    std::map<string,string>::iterator findAttribute(const string &attribute){ return _attributes.find(attribute); }
    /**
     * \brief returns and iterator to an attribute
     */    
    std::map<string,string>::iterator lastAttribute(){ return _attributes.end(); }   
    /**
     * \brief stores output format
     */    
    static int outputFormat(){return _format;}; 
    /**
     * \brief stores output level
     */    
    static int outputLevel(){return _output_level;}; 
    /**
     * \brief stores output indentation
     */        
    static int outputIndent(){return _output_indent;}; 
     
private:        
    map<string,Property*> _map;
    map<string,string> _attributes;
    list<Property> _properties;
    
    string _name;
    string _value;
    string _path;

    static const int _format; 
    static const int _output_level;    
    static const int _output_indent;    

};
  

inline Property &Property::set(const string &key, const string &value)
{
    Property &p = get(key);
    p.value() = value;
    return p;
}

inline Property &Property::add(const string &key, const string &value)
{
    string path = _path;
    if(path != "") path = path + ".";
    _properties.push_back(Property(key, value, path + _name ));
    _map[key] = &(_properties.back());
    return _properties.back();
}

inline bool Property::exists(const string &key)
{
    try { get(key); }
    catch(std::exception &err) { return false;}
    return true;
}
    
bool load_property_from_xml(Property &p, string file);

// TO DO: write a better function for this!!!!
template<>
inline bool Property::as<bool>() const
{
    if(_value == "true" || _value == "TRUE" || _value == "1") return true;
    else return false;
}

template<typename T>
inline T Property::as() const
{
    return lexical_cast<T>(_value, "wrong type in " + _path + "."  + _name + "\n");
}

template<>
inline std::string Property::as<std::string>() const
{
    string tmp(_value);
    boost::trim(tmp);
    return tmp;
}

template<>
inline vec Property::as<vec>() const {
    vector<double> tmp;
    Tokenizer tok(as<string > (), " ,");
    tok.ConvertToVector<double>(tmp);
    if (tmp.size() != 3)
        throw runtime_error("Vector has " + boost::lexical_cast<string > (tmp.size()) + " instead of three entries");
    return vec(tmp[0], tmp[1], tmp[2]);
}

template<>
inline vector<unsigned int> Property::as<vector <unsigned int> >() const {
    vector<unsigned int> tmp;
    Tokenizer tok(as<string > (), " ,");
    tok.ConvertToVector<unsigned int>(tmp);
    return tmp;
}

template<>
inline vector<int> Property::as<vector <int> >() const {
    vector<int> tmp;
    Tokenizer tok(as<string > (), " ,\n\t");
    tok.ConvertToVector<int>(tmp);
    return tmp;
}

template<>
inline vector<double> Property::as<vector <double> >() const {
    vector<double> tmp;
    Tokenizer tok(as<string > (), " ,\n\t");
    tok.ConvertToVector<double>(tmp);
    return tmp;
}

inline bool Property::hasAttribute(const string &attribute) {
    std::map<string,string>::iterator it;
    it = _attributes.find(attribute);
    if ( it == _attributes.end() ) return false;
    return true;
}

template<typename T>
inline T Property::getAttribute(std::map<string,string>::iterator it)
{
    if (it != _attributes.end()) {
        return lexical_cast<T>((*it).second);
    } else {
        throw std::runtime_error("attribute " + (*it).first + " not found\n");
    }
}

template<typename T>
inline T Property::getAttribute(const string &attribute)
{
    std::map<string,string>::iterator it;
    
    it = _attributes.find(attribute);
    
    if (it != _attributes.end()) {
        return lexical_cast<T>(_attributes[attribute], "wrong type in attribute " + attribute + " of element " + _path + "."  + _name + "\n");
    } else {
        throw std::runtime_error("attribute " + attribute + " not found\n");
    }
}
template<typename T>
inline void Property::setAttribute(const string &attribute, const T &value)
{
     _attributes[attribute] = lexical_cast<string>(value, "wrong type to set attribute");
}

inline void throwRuntimeError(string message) {
    
}

}}

#endif	/* _VOTCA_TOOLS_PROPERTY_H */
