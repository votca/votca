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

#include <fstream>
#include <vector>
#include <votca/tools/tokenizer.h>
#include <votca/tools/table.h>
#include <stdexcept>
#include <iostream>
#include <boost/algorithm/string/replace.hpp>
#include <votca/tools/lexical_cast.h>

namespace votca { namespace tools {

using namespace boost;
using namespace std;

void Table::resize(int N, bool preserve)
{
    _x.resize(N, preserve);
    _y.resize(N, preserve);
    _flags.resize(N, preserve);
    if (_has_yerr) {
        _yerr.resize(N, preserve);
    }
}

void Table::Load(string filename)
{
    ifstream in;
    in.open(filename.c_str());
    if(!in)
        throw runtime_error(string("error, cannot open file ") + filename);

    setErrorDetails("file " + filename);
    in >> *this;
    
    in.close();   
}

void Table::Save(string filename) const
{
    ofstream out;
    out.open(filename.c_str());
    if(!out)
        throw runtime_error(string("error, cannot open file ") + filename);

    if (_has_comment) {
        string str = "# " + _comment_line;
        boost::replace_all(str, "\n", "\n# ");
        boost::replace_all(str, "\\n", "\n# ");
        out << str <<endl;
    }
    out << (*this);
    
    out.close(); 
}

void Table::clear(void)
{
    _x.clear();
    _y.clear();
    _flags.clear();
    _yerr.clear();
}

// TODO: this functon is weired, reading occours twice, cleanup!!
// TODO: modify function to work properly, when _has_yerr is true
inline istream &operator>>(istream &in, Table& t)
{
    size_t N;
    bool bHasN=false;
    string line;
    int line_number=0;
   t.clear();
    
    // read till the first data line
    while(getline(in, line)) {
        line_number++;
        string conversion_error =  t.getErrorDetails() + ", line " + boost::lexical_cast<string>(line_number);

        // remove comments and xmgrace stuff
        line = line.substr(0, line.find("#"));
        line = line.substr(0, line.find("@"));
    
        // tokenize string and put it to vector
        Tokenizer tok(line, " \t");
        vector<string> tokens;
        tok.ToVector(tokens);
        
        // skip empty lines
        if(tokens.size()==0) continue;
        
        // if first line is only 1 token, it's the size
        if(tokens.size() == 1) {
            N = lexical_cast<int>(tokens[0], conversion_error);
            bHasN = true;
        }
        // it's the first data line with 2 or 3 entries
        else if(tokens.size() == 2) {
            t.push_back(lexical_cast<double>(tokens[0], conversion_error), lexical_cast<double>(tokens[1], conversion_error), 'i');
        }
        else if(tokens.size() > 2) {
           char flag='i';
           string sflag = tokens.back();
            if(sflag == "i" || sflag == "o" || sflag == "u")
                flag = sflag.c_str()[0];
            t.push_back(lexical_cast<double>(tokens[0], conversion_error), lexical_cast<double>(tokens[1], conversion_error), flag);
        }
        else throw runtime_error("error, wrong table format");                                
    }
    
    // read the rest
    while(getline(in, line)) {
        line_number++;
        string conversion_error =  t.getErrorDetails() + ", line " + boost::lexical_cast<string>(line_number);

        // remove comments and xmgrace stuff
        line = line.substr(0, line.find("#"));
        line = line.substr(0, line.find("@"));
    
        // tokenize string and put it to vector
        Tokenizer tok(line, " \t");
        vector<string> tokens;
        tok.ToVector(tokens);
        
        // skip empty lines
        if(tokens.size()==0) continue;
                    
        // it's a data line
        if(tokens.size() == 2) {            
            t.push_back(lexical_cast<double>(tokens[0], conversion_error), lexical_cast<double>(tokens[1], conversion_error), 'i');
        }
        else if(tokens.size() > 2) {
            char flag='i';
            if(tokens[2] == "i" || tokens[2] == "o" || tokens[2] == "u")
                flag = tokens[2].c_str()[0];
            t.push_back(lexical_cast<double>(tokens[0], conversion_error), lexical_cast<double>(tokens[1], conversion_error), flag);
        }
        // otherwise error
        else throw runtime_error("error, wrong table format");                                
        
        // was size given and did we read N values?
        if(bHasN)
            if(--N == 0) break;
    }

    return in;
}

void Table::GenerateGridSpacing(double min, double max, double spacing)
{
    int n = floor((max - min)/spacing + 1.000000001);
    resize(n);
    int i=0;
    for(double x=min; i<n; x+=spacing, ++i)
        _x[i] = x;
}

void Table::Smooth(int Nsmooth)
{
    while(Nsmooth-- > 0)
        for(int i=1; i<size()-1; ++i)
            _y[i] =0.25*(_y[i-1] + 2*_y[i] +  _y[i+1]);
}
    
}}
