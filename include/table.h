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


// the entry is invalid, e.g. could not be calculated (ln(0), ...)
#define TBL_INVALID    1

/**
    \brief class to store tables like rdfs, tabulated potentials, etc
 
    \think about weather to make this a template, can be used in histogram
    as well, of for counting with integeers...
 
 */
class Table
{
public:       
    Table() {};
    Table(Table &tbl);
    
    ~Table() {};
    
    void resize(int N) { _x.resize(N); _y.resize(N); _flags.resize(N); }
    int size() const {return _x.size(); }

    double &x(int i) { return _x[i]; }
    double &y(int i) { return _y[i]; }
    unsigned short &flags(int i) { return _flags[i]; }

    void set(const int &i, const double &x, const double &y) { _x[i] = x; _y[i]=y; }
    void set(const int &i, const double &x, const double &y, const int &flags) { _x[i] = x; _y[i]=y; _flags[i] = flags; }

    void Load(string filename);
    void Save(string filename) const;       
    
    ub::vector<double> &x() { return _x; }
    ub::vector<double> &y() { return _x; }
    ub::vector<double> &flags() { return _x; }

private:
    ub::vector<double> _x;
    ub::vector<double> _y;       
    ub::vector<unsigned short>   _flags;

    friend ostream &operator<<(ostream &out, const Table& v);
    friend istream &operator>>(istream &out, Table& v);

};

Table::Table(Table &tbl)
{
    resize(tbl.size());
    _x = tbl._x;
    _y = tbl._y;
    _flags = tbl._flags;
}

inline ostream &operator<<(ostream &out, const Table& t)
{
    out << t.size() << endl;
    for(int i=0; i<t._x.size(); ++i) {
        out << t._x[i] << " " << t._y[i] << " " << t._flags[i] << endl;
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
        in >> t._flags[i];
    }
    return in;
}

#endif	/* _TABLE_H */

