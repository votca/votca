/* 
 * File:   average.h
 * Author: ruehle
 *
 * Created on July 17, 2008, 2:08 PM
 */

#ifndef _AVERAGE_H
#define	_AVERAGE_H

// do not use this calss yet!
template<typename T>
class Average
{
public:
    Average();
    ~Average() {}
    
    void Process(T &v);
    template<typename iterator_type>
    void ProcessRange(const iterator_type &begin, const iterator_type   &end);
    double CalcDev();
    double GetAv();
    
private:
    
    T _av; // average
    T _m2; // second moment
    size_t _n;
};

template<typename T>
Average<T>::Average()
: _n(0) {}

template <>
inline Average<double>::Average()
: _n(0), _av(0), _m2(0) {}

template<typename T>
inline void Average<T>::Process(T &value)
{
    if(_n==0){_av = value; _m2= value*value;}
    else {
        _av = _av*(double)_n/(double)(_n+1) + value / (double)(_n+1);
        _n++;
        _m2 += value*value;
    }
}

template<>
inline void Average<double>::Process(double &value)
{
    _av = _av*(double)_n/(double)(_n+1) + value / (double)(_n+1);
    _n++;
    _m2 += value*value;
}

template<typename T>
template<typename iterator_type>
void Average<T>::ProcessRange(const iterator_type &begin, const iterator_type   &end){ 
    for(iterator_type iter=begin; iter!=end; ++iter){
        Process(*iter);
    }
}

template<typename T>
double Average<T>::CalcDev(){
    double dev = 0.0;
    dev = sqrt((_m2-_n*_av*_av)/(_n-1));
    return dev;
}

template<typename T>
double Average<T>::GetAv(){
    return _av;
}
#endif	/* _AVERAGE_H */

