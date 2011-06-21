/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _TABLE_H
#define	_TABLE_H

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <string>

namespace votca { namespace tools {

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
    Table() ;
    Table(Table &tbl);
        
    ~Table() {};
    
    void clear(void);
    
    void GenerateGridSpacing(double min, double max, double spacing);
//    void resize(int N, bool preserve=true) { _x.resize(N, preserve); _y.resize(N, preserve); _flags.resize(N, preserve); }
    void resize(int N, bool preserve=true);
    int size() const {return _x.size(); }

    double &x(int i) { return _x[i]; }
    double &y(int i) { return _y[i]; }
    char &flags(int i) { return _flags[i]; }
    double &yerr(int i) { return _yerr[i]; }

    void set(const int &i, const double &x, const double &y) { _x[i] = x; _y[i]=y; }
    void set(const int &i, const double &x, const double &y, const char &flags) { _x[i] = x; _y[i]=y; _flags[i] = flags; }
    void set(const int &i, const double &x, const double &y, const char &flags, const double &yerr) { _x[i] = x;
                                                                                            _y[i]=y; _flags[i] = flags; _yerr[i] = yerr; }

    void set_comment(const string comment) {_has_comment=true; _comment_line = comment;}

    void Load(string filename);
    void Save(string filename) const;       
    
    void Smooth(int Nsmooth);

    bool GetHasYErr() { return _has_yerr; }
    void SetHasYErr(bool has_yerr) { _has_yerr = has_yerr; }
    
    ub::vector<double> &x() { return _x; }
    ub::vector<double> &y() { return _y; }
    ub::vector<char> &flags() { return _flags; }
    ub::vector<double> &yerr() { return _yerr; }
    
    void push_back(double x, double y, char flags);

    const string &getErrorDetails() {
        return _error_details;
    }

    void setErrorDetails(string str) {
        _error_details = str;
    }

private:
    ub::vector<double> _x;
    ub::vector<double> _y;       
    ub::vector<char>   _flags;
    ub::vector<double> _yerr;
    string _error_details;

    bool _has_yerr;
    bool _has_comment;

    friend ostream &operator<<(ostream &out, const Table& v);
    friend istream &operator>>(istream &out, Table& v);

    string _comment_line;

};

inline Table::Table()
{
    _has_yerr = false;
    _error_details = "";
}

inline Table::Table(Table &tbl)
{
    resize(tbl.size());
    _x = tbl._x;
    _y = tbl._y;
    _flags = tbl._flags;
    _has_yerr = tbl._has_yerr;
    if (_has_yerr) _yerr = tbl._yerr;
    _has_comment = false;
    _error_details = "";
}

inline ostream &operator<<(ostream &out, const Table& t)
{
    if ( t._has_yerr ) {
        for(size_t i=0; i<t._x.size(); ++i) {
            out << t._x[i] << " " << t._y[i] << " " << t._yerr[i] << " " << t._flags[i] << endl;
        }
    }
    else {
        //out << t.size() << endl;
        for(size_t i=0; i<t._x.size(); ++i) {
            out << t._x[i] << " " << t._y[i] << " " << t._flags[i] << endl;
        }
    }
    return out;
}
// TODO: modify this function to be able to treat _has_yerr == true
inline void Table::push_back(double x, double y, char flags)
{
    size_t n=size();
    resize(n+1);
    _x[n] = x;
    _y[n] = y;
    _flags[n] = flags;
}

istream &operator>>(istream &in, Table& t);

}}

#endif	/* _TABLE_H */

