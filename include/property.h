/* 
 * File:   property.h
 * Author: ruehle
 *
 * Created on July 17, 2008, 2:16 PM
 */

#ifndef _PROPERTY_H
#define	_PROPERTY_H

#include <map>
#include <string>
#include <boost/lexical_cast.hpp>
#include <list>
#include <boost/algorithm/string/trim.hpp>
#include <stdexcept>

using namespace std;

/**
    \brief class to manage properties and options, also used as an xml-wrapper

    This class stores names and values in a hierachical structure like it's
    used in an xml format. The structure can be eiter filled manually
    or read in from an xml file using load_property_from_xml.

 */
class Property {
public:
    Property() : _path("") {}
    
    Property(const string &name, const string &value, const string &path) 
        : _name(name), _value(value), _path(path) {}
    
    /// \brief add a new property to structure
    Property &add(const string &key, const string &value);
    /// \brief set value of existing property
    Property &set(const string &key, const string &value);
    
    /// \brief get existing property
    ///
    /// This function tries to find a property specified by key separated
    /// by "." to step down hierarchie. If the property is not
    /// found a runtime_exception is thrown.
    Property &get(const string &key);

    /// \brief check weather property exists
    bool exists(const string &key);
    
    /// \brief select property based on a filter
    ///
    /// returns a list of properties that match the key kriteria including 
    /// wildcards "*" and "?". Example: "base.item*.value"
    ///
    std::list<Property *> Select(const string &filter);
    
    /// \brief reference to value of property
    string &value() { return _value; }
    /// \brief name of property
    string name() { return _name; }
    /// \brief full path of property
    string path() { return _path; }

    /// \brief return value as type
    ///
    /// returns the value after type conversion, e.g.
    /// p.as<int>() returns an integer

    template<typename T>
    T as() const;

    /// \brief does the property has childs?
    bool HasChilds() { return !_map.empty(); }

    /// iterator to iterate over properties
    typedef list<Property>::iterator iterator;
    
    /// \brief iterator to first child property
    iterator begin() { return _properties.begin(); }
    /// \brief end iterator for child properties
    iterator end() { return _properties.end(); }
    /// \brief number of child properties
    list<Property>::size_type size() { return _properties.size(); }

private:        
    map<string,Property*> _map;
    list<Property> _properties;
    
    string _name;
    string _value;
    string _path;
    
    static void PrintNode(std::ostream &out, const string &prefix, Property &p);
    
    friend std::ostream &operator<<(std::ostream &out, Property& p);
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
    catch(std::exception err) { return false;}
    return true;
}
    
bool load_property_from_xml(Property &p, string file);
std::ostream &operator<<(std::ostream &out, Property& p);

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
    try {
       return boost::lexical_cast<T>(_value);
    }
    catch(boost::bad_lexical_cast error ) {
        throw std::runtime_error("wrong type in " + _path + "."  + _name + "\n" + error.what());
    }
}

template<>
inline std::string Property::as<std::string>() const
{
    string tmp(_value);
    boost::trim(tmp);
    return tmp;
}

#endif	/* _PROPERTY_H */
