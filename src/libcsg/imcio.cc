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

#include <boost/algorithm/string/trim.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <votca/tools/rangeparser.h>
#include <votca/tools/table.h>
#include <boost/lexical_cast.hpp>
#include <votca/tools/tokenizer.h>
#include <votca/tools/getline.h>
#include <iostream>
#include <votca/csg/imcio.h>

namespace votca { namespace csg {

typedef Eigen::MatrixXd group_matrix;
using namespace std;

void imcio_write_dS(const string &file, Eigen::VectorXd &r, Eigen::VectorXd &dS, std::list<int> *list)
{
    // write the dS
    ofstream out_dS;
    out_dS.open(file.c_str());
    out_dS << setprecision(8);
    if (!out_dS)
        throw runtime_error(string("error, cannot open file ") + file);

    if(list == NULL) {
        for (int i = 0; i < dS.size(); ++i) {
            out_dS << r[i] << " " << dS[i] << endl;
        }
    }
    else {
        for(std::list<int>::iterator i = list->begin(); i!=list->end(); ++i) {
            out_dS << r[*i] << " " << dS[*i] << endl;
        }
    }


    out_dS.close();
    cout << "written " << file << endl;
}

void imcio_write_matrix(const string &file, Eigen::MatrixXd &gmc, std::list<int> *list)
{
    ofstream out_A;
    out_A.open(file.c_str());
    out_A << setprecision(8);

    if (!out_A)
        throw runtime_error(string("error, cannot open file ") + file);

    if(list == NULL) {
        for (int i = 0; i < gmc.rows(); ++i) {
            for (int j = 0; j < gmc.cols(); ++j) {
                out_A << gmc(min(i,j), max(i,j)) << " ";
            }
            out_A << endl;
        }
    }
    else {
        for(std::list<int>::iterator i = list->begin(); i!=list->end(); ++i) {
            for(std::list<int>::iterator j = list->begin(); j!=list->end(); ++j) {
                out_A << gmc(*i, *j) << " ";
            }
            out_A << endl;
        }
    }
out_A.close();
    cout << "written " << file << endl;
}

void imcio_write_index(const string &file, vector<string> &names, vector<RangeParser> &ranges)
{
    // write the index

    ofstream out_idx;
    out_idx.open(file.c_str());

    if (!out_idx)
        throw runtime_error(string("error, cannot open file ") + file);
    
    for (size_t i = 0; i < names.size(); ++i)
        out_idx << names[i] << " " << ranges[i] << endl;
    
    out_idx.close();
    cout << "written " << file << endl;
}

void imcio_read_dS(const string &filename, Eigen::VectorXd &r, Eigen::VectorXd &dS)
{
    Table tbl;
    tbl.Load(filename);
    
    r.resize(tbl.size());
    dS.resize(tbl.size());
    
    for( int i=0; i<tbl.size(); ++i) {
        r(i) = tbl.x(i);
        dS(i) = tbl.y(i);
    }
}

void imcio_read_matrix(const string &filename, Eigen::MatrixXd &gmc)
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
        if(!is_initialized)
            gmc.resize(tokens.size(),tokens.size());
        is_initialized=true;

        if(gmc.rows()!=int(tokens.size()))
            throw runtime_error(string("error loading ")
                    + filename + ": size mismatchm, number of columns differ");
        for(size_t i=0; i<tokens.size(); ++i)
            gmc(line_count,i) = stod(tokens[i]);
        ++line_count;
    }
    if(line_count != gmc.rows())
            throw runtime_error(string("error loading ")
                    + filename + ": size mismatch, not enough lines");
    in.close();
}

void imcio_read_index(const string &filename, vector<string> &names, vector<RangeParser> &ranges)
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
        
        RangeParser rp;
        rp.Parse(range);
        names.push_back(name);
        ranges.push_back(rp);

    }
    in.close();
}

}}
