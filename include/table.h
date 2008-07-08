/* 
 * File:   table.h
 * Author: victorr
 *
 * Created on June 11, 2008, 1:34 PM
 */

#ifndef _TABLE_H
#define	_TABLE_H

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
using namespace std;

namespace ub = boost::numeric::ublas;

/**
    \brief class to store tables like rdfs, tabulated potentials, etc
 
    \think about weather to make this a template, can be used in histogram
    as well, of for counting with integeers...
 
 */
class Table
{
public:       
    Table() {};
    ~Table() {};
    
    void resize(int N) { _x.resize(N); _y.resize(N); }
    int size() const {return _x.size(); }
    
    ub::vector<double> &x() { return _x; }
    ub::vector<double> &y() { return _y; }
    double &x(int i) { return _x[i]; }
    double &y(int i) { return _y[i]; }
    
    void set(const int &i, const double &x, const double &y) { _x[i] = x; _y[i]=y; }

    void Load(string filename);
    void Save(string filename) const;       
    
    void ScaleRdf();
    
private:
    ub::vector<double> _x;
    ub::vector<double> _y;       
    
    friend ostream &operator<<(ostream &out, const Table& v);
    friend istream &operator>>(istream &out, Table& v);

};

inline ostream &operator<<(ostream &out, const Table& t)
{
    out << t.size() << endl;
    for(int i=0; i<t._x.size(); ++i) {
        out << t._x[i] << " " << t._y[i] << endl;
    }
    return out;
}

inline istream &operator>>(istream &in, Table& t)
{
    int N;
    in >> N;
    t.resize(N);
    for(int i=0; i<N; ++i) {
        in >> t._x[i];
        in >> t._y[i];
    }
    return in;
}

#endif	/* _TABLE_H */

