/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_FILEFORMATFACTORY_H
#define	_VOTCA_CSG_FILEFORMATFACTORY_H

#include <string>
#include <iostream>
#include <votca/tools/objectfactory.h>
#include <votca/tools/tokenizer.h>

namespace votca { namespace csg {
using namespace votca::tools;

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

}}

#endif	/* _VOTCA_CSG_FILEFORMATFACTORY_H */

