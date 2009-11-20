/* 
 * File:   table.cc
 * Author: victorr
 *
 * Created on June 11, 2008, 1:48 PM
 */

#include <fstream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "tokenizer.h"
#include "table.h"
#include <stdexcept>

using namespace boost;

void Table::Load(string filename)
{
    ifstream in;
    in.open(filename.c_str());
    if(!in)
        throw runtime_error(string("error, cannot open file ") + filename);
    
    in >> *this;
    
    in.close();   
}

void Table::Save(string filename) const
{
    ofstream out;
    out.open(filename.c_str());
    if(!out)
        throw runtime_error(string("error, cannot open file ") + filename);
    
    out << (*this);
    
    out.close(); 
}

void Table::clear(void)
{
    _x.clear();
    _y.clear();
    _flags.clear();
}

// TODO: this functon is weired, reading occours twice, cleanup!!
inline istream &operator>>(istream &in, Table& t)
{
    size_t N;
    bool bHasN=false;
    string line;
    
   t.clear();
    
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
        
        // if first line is only 1 token, it's the size
        if(tokens.size() == 1) {
            N = lexical_cast<int>(tokens[0]);
            bHasN = true;
        }
        // it's the first data line with 2 or 3 entries
        else if(tokens.size() == 2) {
            t.push_back(lexical_cast<double>(tokens[0]), lexical_cast<double>(tokens[1]), 'i');            
        }
        else if(tokens.size() > 2) {
           char flag='i';
            if(tokens[2] == "i" || tokens[2] == "o" || tokens[2] == "u")
                flag = tokens[2].c_str()[0];
            t.push_back(lexical_cast<double>(tokens[0]), lexical_cast<double>(tokens[1]), flag);
        }
        else throw runtime_error("error, wrong table format");                                
    }
    
    // read the rest
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
                    
        // it's a data line
        if(tokens.size() == 2) {            
            t.push_back(lexical_cast<double>(tokens[0]), lexical_cast<double>(tokens[1]), 'i');            
        }
        else if(tokens.size() > 2) {
            char flag='i';
            if(tokens[2] == "i" || tokens[2] == "o" || tokens[2] == "u")
                flag = tokens[2].c_str()[0];
            t.push_back(lexical_cast<double>(tokens[0]), lexical_cast<double>(tokens[1]), flag);
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
    
