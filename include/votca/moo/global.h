#ifndef FILE_GLOB
#define FILE_GLOB

#include <iostream>

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

#endif //FILE_GLOB
