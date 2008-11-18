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
#include <tokenizer.h>

using namespace std;

class Property {
public:
    Property &set(const string &key, const string &value);
    Property &get(const string &key);
    
    string &value() { return _value; }
    template<typename T>
    T &as() { return boost::lexical_cast<T>(_value); }

    bool HasChilds() { return !_map.empty(); }

    typedef map<string, Property>::iterator iterator;
    iterator begin() { return _map.begin(); }
    iterator end() { return _map.end(); }
    map<string,Property>::size_type size() { return _map.size(); }

private:
    map<string,Property> _map;
    string _value;
    
    static void PrintNode(std::ostream &out, const string &prefix, Property &p);
    
    friend std::ostream &operator<<(std::ostream &out, Property& p);
};

inline Property &Property::set(const string &key, const string &value)
{
    _map[key].value() = value;
    return _map[key];
}

bool load_property_from_xml(Property &p, string file);
std::ostream &operator<<(std::ostream &out, Property& p);

#endif	/* _PROPERTY_H */
