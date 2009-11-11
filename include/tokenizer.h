// 
// File:   tools.h
// Author: ruehle
//
// Created on May 11, 2007, 12:37 PM
//

#ifndef _tools_H
#define	_tools_H

#include <string>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class Tokenizer 
{                        
public:
    typedef boost::tokenizer<boost::char_separator<char> >::iterator
        iterator;
    
    Tokenizer(const std::string &str, const char *separators) {
        boost::char_separator<char> sep(separators);
        tok = new boost::tokenizer<boost::char_separator<char> >(str, sep);
        //boost::escaped_list_separator<char> sep(" ", separators, "\"");
        //tok = new boost::tokenizer<boost::escaped_list_separator<char> >(str, sep);
        
    }
    
    ~Tokenizer() {
        delete tok;
    }
    
    iterator begin() { return tok->begin(); }
    iterator end() { return tok->end(); }
    
    void ToVector(std::vector<std::string> &v) {
        for(iterator iter=begin(); iter!=end(); ++iter)
            v.push_back(*iter);
    }

    template < typename T >
    void ConvertToVector(std::vector<T> &v){
        std::vector<std::string> tmp;
        ToVector(tmp);
        v.resize(tmp.size());
        transform(tmp.begin(), tmp.end(), v.begin(),
            boost::lexical_cast<T, std::string>);
    }
    
private:
    boost::tokenizer< boost::char_separator<char> > *tok;
};


int wildcmp(const char *wild, const char *string);

#endif	/* _tools_H */

