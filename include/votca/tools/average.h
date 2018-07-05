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

#ifndef _AVERAGE_H
#define	_AVERAGE_H

namespace votca { namespace tools {

// do not use this calss yet!
template<typename T>
class Average
{
public:
    Average();
    ~Average() {}
    
    void Process(const T &v);
    void Clear();
    template<typename iterator_type>
    void ProcessRange(const iterator_type &begin, const iterator_type   &end);
    
    T CalcDev();
    T CalcSig2();
    const T &getAvg();
    const T getM2();
    size_t getN();
    
private:
    
    size_t _n;
    T _av; // average
    T _m2; // second moment
};

template<typename T>
Average<T>::Average()
: _n(0) {}

template <>
inline Average<double>::Average()
: _n(0), _av(0), _m2(0) {}

template<typename T>
inline void Average<T>::Process(const T &value)
{ 
   _av = _av*(double)_n/(double)(_n+1) + value / (double)(_n+1);
   _n++;
   _m2 += value*value;   
}

template<typename T>
inline void Average<T>::Clear()
{
   _av = 0;
   _n = 0;
   _m2 = 0 ;
}

/*
template<>
inline void Average<double>::Process(const double &value)
{
    _av = _av*(double)_n/(double)(_n+1) + value / (double)(_n+1);
    _n++;
    _m2 += value*value;
}
*/
template<typename T>
template<typename iterator_type>
void Average<T>::ProcessRange(const iterator_type &begin, const iterator_type   &end){ 
    for(iterator_type iter=begin; iter!=end; ++iter){
        Process(*iter);
    }
}

template<typename T>
T Average<T>::CalcDev(){
    double dev = 0.0;
    dev = sqrt((_m2-_n*_av*_av)/(_n-1));
    return dev;
}

template<typename T>
T Average<T>::CalcSig2(){
    double dev = 0.0;
    dev = _m2/_n-_av*_av ;
    return dev;
}

template<typename T>
const T &Average<T>::getAvg(){
    return _av;
}

template<typename T>
const T Average<T>::getM2(){
    double m2 = 0.0;
    m2 = _m2/_n;
    return m2;
}

template<typename T>
size_t Average<T>::getN(){
    return _n;
}

}}

#endif	/* _AVERAGE_H */

