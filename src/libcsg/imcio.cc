/* 
 * File:   imcio.cpp
 * Author: ruehle
 *W
 * Created on September 14, 2009, 5:36 PM
 */

#include <boost/algorithm/string/trim.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <tools/rangeparser.h>
#include <tools/table.h>
#include <boost/lexical_cast.hpp>
#include <tools/tokenizer.h>
#include <iostream>
#include "imcio.h"

typedef ub::symmetric_matrix<double> group_matrix;
using namespace std;

void imcio_write_dS(string file, ub::vector<double> &r, ub::vector<double> &dS)
{
    // write the dS
    ofstream out_dS;
    out_dS.open(file.c_str());
    out_dS << setprecision(8);
    if (!out_dS)
        throw runtime_error(string("error, cannot open file ") + file);

    for (int i = 0; i < dS.size(); ++i) {
        out_dS << r[i] << " " << dS[i] << endl;
    }

    out_dS.close();
    cout << "written " << file << endl;
}

void imcio_write_matrix(string file, ub::symmetric_matrix<double> gmc)
{
    ofstream out_A;
    out_A.open(file.c_str());
    out_A << setprecision(8);

    if (!out_A)
        throw runtime_error(string("error, cannot open file ") + file);

    for (group_matrix::size_type i = 0; i < gmc.size1(); ++i) {
        for (group_matrix::size_type j = 0; j < gmc.size2(); ++j) {
            out_A << gmc(i, j) << " ";
        }
        out_A << endl;
    }
    out_A.close();
    cout << "written " << file << endl;
}

void imcio_write_index(string file, vector<string> &names, vector<RangeParser> &ranges)
{
    // write the index

    ofstream out_idx;
    out_idx.open(file.c_str());

    if (!out_idx)
        throw runtime_error(string("error, cannot open file ") + file);
    
    for (int i = 0; i < names.size(); ++i)
        out_idx << names[i] << " " << ranges[i] << endl;
    
    out_idx.close();
    cout << "written " << file << endl;
}

void imcio_read_dS(string filename, ub::vector<double> &r, ub::vector<double> &dS)
{
    Table tbl;
    tbl.Load(filename);
    
    r.resize(tbl.size());
    dS.resize(tbl.size());
    
    for(int i=0; i<tbl.size(); ++i) {
        r(i) = tbl.x(i);
        dS(i) = tbl.y(i);
    }
}

void imcio_read_matrix(string filename, ub::symmetric_matrix<double> gmc)
{
    ifstream in;
    in.open(filename.c_str());

    bool is_initialized = false;
    if(!in)
        throw runtime_error(string("error, cannot open file ") + filename);

    int line_count =0;
    string line;
    // read till the first data line
    while(getline(in, line)) {
        // remove comments and xmgrace stuff
        line = line.substr(0, line.find("#"));
        line = line.substr(0, line.find("@"));

        // tokenize string and put it to vector
        Tokenizer tok(line, " \t");
        vector<string> tokens;
        tok.ToVector(tokens);

        // skip empty lines
        if(tokens.size()==0) continue;

        gmc.resize(tokens.size());
        is_initialized=true;

        if(gmc.size1()!=tokens.size())
            throw runtime_error(string("error loading ")
                    + filename + ": size mismatchm, number of columns differ");

        for(int i=0; i<tokens.size(); ++i)
            gmc(line_count,i) = boost::lexical_cast<double>(tokens[i]);
        ++line_count;
    }
    if(line_count != gmc.size1())
            throw runtime_error(string("error loading ")
                    + filename + ": size mismatch, not enough lines");
    in.close();
}

void imcio_read_index(string filename, vector<string> &names, vector<RangeParser> &ranges)
{
    ifstream in;
    in.open(filename.c_str());
    if(!in)
        throw runtime_error(string("error, cannot open file ") + filename);

    names.clear();
    ranges.clear();
    string line;
    
    // read till the first data line
    while(getline(in, line)) {
        // remove comments and xmgrace stuff
        line = line.substr(0, line.find("#"));
        line = line.substr(0, line.find("@"));

        boost::trim(line);

        size_t found;
        found = line.find(" ");
        if(found==string::npos)
            throw runtime_error(string("wrong format in ") + filename);

        string name = line.substr(0, found);

        string range = line.substr(found);
        cout << "<" << name << "><" << range << ">\n";

        RangeParser rp;
        rp.Parse(range);
        names.push_back(name);
        ranges.push_back(rp);

    }
    in.close();
}