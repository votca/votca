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
    
    void clear(void);
    
    void GenerateGridSpacing(double min, double max, double spacing);
    void resize(int N, bool preserve=true) { _x.resize(N, preserve); _y.resize(N, preserve); _flags.resize(N, preserve); }
    int size() const {return _x.size(); }

    double &x(int i) { return _x[i]; }
    double &y(int i) { return _y[i]; }
    char &flags(int i) { return _flags[i]; }

    void set(const int &i, const double &x, const double &y) { _x[i] = x; _y[i]=y; }
    void set(const int &i, const double &x, const double &y, const char &flags) { _x[i] = x; _y[i]=y; _flags[i] = flags; }

    void Load(string filename);
    void Save(string filename) const;       
    
    void Smooth(int Nsmooth);
    
    ub::vector<double> &x() { return _x; }
    ub::vector<double> &y() { return _y; }
    ub::vector<char> &flags() { return _flags; }
    
    void push_back(double x, double y, char flags);

private:
    ub::vector<double> _x;
    ub::vector<double> _y;       
    ub::vector<char>   _flags;

    friend ostream &operator<<(ostream &out, const Table& v);
    friend istream &operator>>(istream &out, Table& v);

};

inline Table::Table(Table &tbl)
{
    resize(tbl.size());
    _x = tbl._x;
    _y = tbl._y;
    _flags = tbl._flags;
}

inline ostream &operator<<(ostream &out, const Table& t)
{
    //out << t.size() << endl;
    for(int i=0; i<t._x.size(); ++i) {
        out << t._x[i] << " " << t._y[i] << " " << t._flags[i] << endl;
    }
    return out;
}

inline void Table::push_back(double x, double y, char flags)
{
    size_t n=size();
    resize(n+1);
    _x[n] = x;
    _y[n] = y;
    _flags[n] = flags;
}

inline istream &operator>>(istream &in, Table& t);

#endif	/* _TABLE_H */

