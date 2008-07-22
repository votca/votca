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
    
private:
    
    T _av;
    size_t n;
};

template<typename T>
void Average::Average()
  : n(0)
{}

template<typename T>
void Average::Process(T &v)
{
    v = _av*(double)n/(double)(n+1) + v / (double)(n+1);
    ++n;
}

#endif	/* _AVERAGE_H */

