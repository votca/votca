// 
// File:   fileformatfactory.h
// Author: ruehle
//
// Created on January 11, 2008, 4:28 PM
//

#ifndef _FILEFORMATFACTORY_H
#define	_FILEFORMATFACTORY_H

#include <string>
#include <iostream>
#include <votca/tools/objectfactory.h>
#include <votca/tools/tokenizer.h>

using namespace std;

template<typename T>
class FileFormatFactory
: public ObjectFactory<std::string, T>
{
public:
    FileFormatFactory() {}
    
    T *Create(const string &file);
};

template<typename T>
T *FileFormatFactory<T>::Create(const string &file)
{
    string filetype = "";
    Tokenizer tok(file, ".");       
    for(Tokenizer::iterator iter=tok.begin();iter!=tok.end();iter++)
        filetype = *iter;
    try {
        return ObjectFactory<string,T>::Create(filetype);
    } catch(std::exception &error) {}
    return NULL;
}

#endif	/* _FILEFORMATFACTORY_H */

