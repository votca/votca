// 
// File:   tabulatedpotential.cc
// Author: ruehle
//
// Created on August 2, 2007, 3:18 PM
//

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <boost/lexical_cast.hpp>
#include "analysistool.h"
#include <votca/tools/histogram.h>
#include "version.h"
#include "bondedstatistics.h"
#include "tabulatedpotential.h"

using namespace std;
using namespace boost;

TabulatedPotential::TabulatedPotential()
{
    _tab_smooth1 = _tab_smooth2 = 0;
    _T = 300;
}

void TabulatedPotential::Register(map<string, AnalysisTool *> &lib)
{
    lib["tab"] = this;
    lib["hist"] = this;
}

void TabulatedPotential::Command(BondedStatistics &bs, string cmd, vector<string> &args)
{
    if(args[0] == "set") {
        if(cmd == "hist") SetOption(_hist_options, args);
        else if(cmd == "tab") {
            if(!SetOption(_tab_options, args)) {
                if(args.size() >2) {
                    if(args[1] == "smooth1")
                        _tab_smooth1 = lexical_cast<int>(args[2]);       
                    else if(args[1] == "smooth2")
                        _tab_smooth2 = lexical_cast<int>(args[2]);       
                    else if(args[1] == "T")
                        _T = lexical_cast<double>(args[2]);
                    else {
                        cout << "unknown option " << args[2] << endl;
                        return;
                    }
                }
            }        
            if(args.size() <=2) {
                cout << "smooth1: " << _tab_smooth1 << endl;
                cout << "smooth2: " << _tab_smooth2 << endl;
                cout << "T: " << _T << endl;
            }
        }
    }
    else if(args.size() >= 2) {
        if(cmd == "hist") WriteHistogram(bs, args);
        else if(cmd == "tab") WritePotential(bs, args);
    }
    else cout << "wrong number of arguments" << endl;
}

bool TabulatedPotential::SetOption(Histogram::options_t &op, const vector<string> &args)
{
    if(args.size() >2) {
        if(args[1] == "n")
            op._n = lexical_cast<int>(args[2]);       
        else if(args[1] == "min") {
            op._min = lexical_cast<double>(args[2]);
        }
        else if(args[1] == "max")
            op._max = lexical_cast<double>(args[2]);
        else if(args[1] == "periodic")
            op._periodic = lexical_cast<bool>(args[2]);
        else if(args[1] == "auto")
            op._auto_interval = lexical_cast<bool>(args[2]);
        else if(args[1] == "extend")
            op._extend_interval = lexical_cast<bool>(args[2]);
        else if(args[1] == "normalize")
            op._normalize = lexical_cast<bool>(args[2]);
        else if(args[1] == "scale") {
            if(args[2]=="no" || args[2]=="bond" || args[2]=="angle")
                op._scale = args[2];
            else {
                cout << "scale can be: no, bond or angle\n";
            }
        }
        else {
            return false;        
        }
    }
    else {
        cout << "n: " << op._n << endl;
        cout << "min: " << op._min << endl;
        cout << "max: " << op._max << endl;
        cout << "periodic: " << op._periodic << endl;
        cout << "auto: " << op._auto_interval << endl;
        cout << "extend: " << op._extend_interval << endl;
        cout << "scale: " << op._scale << endl;
        cout << "normalize: " << op._normalize << endl;
    }
    return true;
}

void TabulatedPotential::WriteHistogram(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = bs.BondedValues().select(args[i], sel);
    Histogram h(_hist_options);
    h.ProcessData(sel);
    out.open(args[0].c_str());
/*    out << "# histogram, created csg version " <<  VERSION_STR  << endl;
    out << "# n = " << _hist_options._n << endl;
    out << "# min = " << _hist_options._min << endl;
    out << "# max = " << _hist_options._max << endl;
    out << "# periodic = " << _hist_options._periodic << endl;
    out << "# auto = " << _hist_options._auto_interval << endl;
    out << "# extend = " << _hist_options._extend_interval << endl;
    out << "# scale = " << _hist_options._scale << endl;*/
    out << h ;
    out.close();
    cout << "histogram created using " << sel->size() << " data-rows, written to " << args[0] << endl;    
    delete sel;
}


void TabulatedPotential::CalcForce(vector<double> &U, vector<double> &F, double dx, bool bPeriodic)
{
    size_t n=U.size();
    double f = 0.5/dx;
    F.resize(n);
    if(bPeriodic)
        F[n-1] = F[0] = -(U[1] - U[n-2])*f;
    else {
        F[0] = -(U[1] - U[0])*2*f;
        F[n-1] = -(U[n-1] - U[n-2])*2*f;
    }
    for(size_t i=1; i<n-1; i++)
        F[i] = -(U[i+1] - U[i-1])*f;
}

void TabulatedPotential::WritePotential(BondedStatistics &bs, vector<string> &args)
{
   ofstream out;
    DataCollection<double>::selection *sel = NULL;

    for(size_t i=1; i<args.size(); i++)
        sel = bs.BondedValues().select(args[i], sel);
    Histogram h(_tab_options);
    h.ProcessData(sel);
    for(int i=0; i<_tab_smooth1; ++i)
        Smooth(h.getPdf(), _tab_options._periodic);
    BoltzmannInvert(h.getPdf(), _T);
    for(int i=0; i<_tab_smooth2; ++i)
        Smooth(h.getPdf(), _tab_options._periodic);
    out.open(args[0].c_str());
    
/*       out << "# tabulated potential, created csg version " VERSION_STR  << endl;
    out << "# n = " << _tab_options._n << endl;
    out << "# min = " << _tab_options._min << endl;
    out << "# max = " << _tab_options._max << endl;
    out << "# periodic = " << _tab_options._periodic << endl;
    out << "# auto = " << _tab_options._auto_interval << endl;
    out << "# extend = " << _tab_options._extend_interval << endl;
    out << "# scale = " << _tab_options._scale << endl;
    out << "# smooth1 = " << _tab_smooth1 << endl;
    out << "# smooth2 = " << _tab_smooth2 << endl;
    out << "# T = " << _T << endl;*/

    vector<double> F;
    
    CalcForce(h.getPdf(), F, h.getInterval(), _tab_options._periodic);
    for(int i=0; i<h.getN(); i++) {
        out << h.getMin() + h.getInterval()*((double)i) << " " << h.getPdf()[i] << " " << F[i] << endl;
    }
    out.close();
    cout << "histogram created using " << sel->size() << " data-rows, written to " << args[0] << endl;
    delete sel;
}

void TabulatedPotential::Smooth(vector<double> &data, bool bPeriodic)
{
    double old[3];
    int n=data.size();
    if(bPeriodic) {
        old[0] = data[n-3];
        old[1] = data[n-2];
    }
    else {
        old[0] = data[0];
        old[1] = data[0];
    }
    size_t i;
    for(i=0; i<data.size()-2; i++) {
        old[3] = data[i];
        data[i] = (old[0] + 2.*old[1] + 3.*data[i] + 2.*data[i+1] + data[i+2])/9.;
        old[0]=old[1];
        old[1]=old[3];;
    }
    if(bPeriodic) {
        data[i] = (old[0] + 2.*old[1] + 3.*data[i] + 2.*data[i+1] + data[0])/9.;
        old[0]=old[1];old[1]=data[i];
        data[n-1] = data[0];
    }
    else {
        data[i] = (old[0] + 2.*old[1] + 3.*data[i] + 3.*data[i+1])/9.;
        old[0]=old[1];old[1]=data[i];
        i++;
        data[i] = (old[0] + 2.*old[1] + 6.*data[i])/9.;    
    }
}

void TabulatedPotential::BoltzmannInvert(vector<double> &data, double T)
{
    double _min, _max;
    
    _min = numeric_limits<double>::max();
    _max = numeric_limits<double>::min();
    
    for(size_t i=0; i<data.size(); i++) {
        _max = max(data[i], _max);
        if(data[i] > 0) _min = min(data[i], _min);
    }
    _max = -8.3109*T*log(_max)*0.001;
    _min = -8.3109*T*log(_min)*0.001-_max;
    
    for(size_t i=0; i<data.size(); i++) {
        if(data[i] == 0) data[i] = _min;
        else
            data[i] = -8.3109*T*log(data[i])*0.001 - _max;
    }
}
