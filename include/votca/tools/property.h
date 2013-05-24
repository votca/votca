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

#ifndef _PROPERTY_H
#define	_PROPERTY_H

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

using namespace std;

enum formats{ formatXML, formatTXT, formatLOG, formatTEX, formatT2T, formatHLP };
    
/**
  * \brief class to manage properties and options, also used as an xml-wrapper
  *
  * This class stores names and values in a hierachical (tree like) structure like it's
  * used in an xml format. The structure can be eiter filled manually
  * or read in from an xml file using load_property_from_xml.
  */
class Property {
    
    /// \brief outputs the property to the ostream
    friend std::ostream &operator<<(std::ostream &out, Property& p);
    /// \brief output the content in the xml format
    friend void PrintNodeXML(std::ostream &out, const string &prefix, Property &p, string offset);
    friend void PrintNodeTXT(std::ostream &out, const string &prefix, Property &p);

   
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
     * \brief set an attribute
     */
    template<typename T>
    void setAttribute(const string &attribute, const T &value);

    /**
     * \brief return true if a node has attributes
     */
    bool hasAttributes() { return _attributes.size() > 0; }
    
    /**
     * \brief stores output format
     */    
    static int GetFormat(){return _format;}; 

    /**
     * \brief copies values of a property
     * @param prefix starting path (name0.name1.name2 ...) 
     */
     void CopyValues(string prefix, Property &p);
    
private:        
    map<string,Property*> _map;
    map<string,string> _attributes;
    list<Property> _properties;
    
    string _name;
    string _value;
    string _path;

    static const int _format;    
   
    friend std::ostream &operator<<(std::ostream &out, const Property& p);
    
    
    struct PropertyStackEntry_t {
        Property *prop;
        string comment;
    };
/*
stack<Property *> -> stack< Propertz_stack_entry_t>
*/
    
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

template<typename T>
inline T Property::getAttribute(const string &attribute)
{
    return lexical_cast<T>(_attributes[attribute], "wrong type in attribute " + attribute + " of element " + _path + "."  + _name + "\n");
}

template<typename T>
inline void Property::setAttribute(const string &attribute, const T &value)
{
     _attributes[attribute] = lexical_cast<string>(value, "wrong type to set attribute");
}

inline void throwRuntimeError(string message) {
    
}


class PropertyFormat {
    
public:
    explicit PropertyFormat(formats fmt) : _fmt(fmt){}    
    friend std::ostream& operator << (std::ostream& os, const PropertyFormat& pf)
    {
        os.iword(Property::GetFormat()) = pf._fmt;
        return os;
    }
private:
    int _fmt;     
};

extern PropertyFormat XML;
extern PropertyFormat TXT;
extern PropertyFormat T2T;
extern PropertyFormat LOG;
extern PropertyFormat TEX;
extern PropertyFormat HLP;

}}

#endif	/* _PROPERTY_H */
