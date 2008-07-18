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
    void set(const string &key, const string &value);
    Property &get(const string &key);
    
    string &value() { return _value; }
    template<typename T>
    T &as() { return boost::lexical_cast<T>(_value); }

private:
    map<string,Property> _map;
    string _value;
};

inline void Property::set(const string &key, const string &value)
{
    _map[key].value() = value;
}

Property &Property::get(const string &key)
{
    Tokenizer tok(key, "//");
    Tokenizer::iterator n;
    
    n = tok.begin();
    
    map<string, Property>::iterator iter;
    iter = _map.find(*n);
    if(iter == _map.end())
        throw string("property not found: ") + key;
    
    Property *p(&((*iter).second));
    ++n;
    try {
        for(; n!=tok.end(); ++n)
        p = &p->get(*n);
    }
    catch(string err) {
        throw string("property not found: ") + key;
    }
    
    return *p;
}

#endif	/* _PROPERTY_H */
