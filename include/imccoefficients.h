/* 
 * File:   imccrosscor.h
 * Author: victorr
 *
 * Created on June 10, 2008, 12:01 PM
 */

#ifndef _IMCCROSSCOR_H
#define	_IMCCROSSCOR_H

#include <vector>
#include <ostream>

using namespace std;

class IMCCoefficients
{
public:
    IMCCoefficients() {};
    ~IMCCoefficients() {};
    
    void setN(int N) { _N=N; _cross.resize(N*(N+1)/2); _dist.resize(N); clear(); };
    int getN() { return _N; }
    
    void clear();
    
    void Process(vector<int> dist);
    
    void Average();
    
    void OutputCross(ostream &out);
    void OutputDist(ostream &out);
    
    
protected:
    vector<double> _cross;
    vector<double> _dist;
    int _N;
    int _c;
};

#endif	/* _IMCCROSSCOR_H */

