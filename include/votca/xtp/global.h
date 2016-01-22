/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef FILE_GLOB
#define FILE_GLOB

#include <iostream>

namespace votca { namespace xtp {

using namespace std;

template<typename T>
inline void clearListList(T&obj)
{
    typename T::iterator it = obj.begin();
    for ( ; it != obj.end() ; ++it){
        it->clear();
    }
    obj.clear();
}

template<typename T>
inline void safe_delete(T& obj)
{
        typename T::iterator it =obj.begin();
        for( ; it != obj.end() ; ++it ){
//                cout << "call safe_delete" <<endl;
                delete *it;
 //               cout << "called safe_delete" <<endl;
        }
  //      cout << "About to clear safedelete" <<endl;
        obj.clear();
  //      cout << "Cleared safedelete" <<endl;
}

}}

#endif //FILE_GLOB
