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

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <boost/lexical_cast.hpp>
#include "analysistool.h"
#include <votca/tools/histogram.h>
#include <votca/csg/version.h>
#include "bondedstatistics.h"
#include "tabulatedpotential.h"

using namespace std;
using namespace boost;

/* Boltzmann constant in units [ kJ/(mol K) ] */
const double kB = 0.0083109;

/******************************************************************************
 * Public Facing Methods
 ******************************************************************************/
TabulatedPotential::TabulatedPotential()
{
    _tab_smooth1 = _tab_smooth2 = 0;
    _Temperature = 300;
}

void TabulatedPotential::Register(map<string, AnalysisTool *> &lib)
{
    lib["tab"] = this;
    lib["hist"] = this;
}

void TabulatedPotential::Command(BondedStatistics &bs, string cmd, vector<string> &args)
{
  if(args[0]=="set") {
    if(cmd=="hist") {
      SetOption_(_hist_options, args);
    }else if(cmd=="tab") {
      if(!SetOption_(_tab_options, args)) {
        if(args.size() >2) {
          if(args[1]=="smooth_pdf"){
            _tab_smooth1 = lexical_cast<int>(args[2]);       
          }else if(args[1]=="smooth_pot"){
            _tab_smooth2 = lexical_cast<int>(args[2]);       
          } else if(args[1]=="T"){
            _Temperature = lexical_cast<double>(args[2]);
          }else {
            cout << "unknown option " << args[2] << endl;
            return;
          }
        }
      }        
      if(args.size() <=2) {
        cout << "smooth_pdf: " << _tab_smooth1 << endl;
        cout << "smooth_pot: " << _tab_smooth2 << endl;
        cout << "T: " << _Temperature << endl;
      }
    }
  }else if(args.size() >= 2) {
    if(cmd=="hist"){ 
      WriteHistogram(bs, args);
    }else if(cmd=="tab"){ 
      WritePotential(bs, args);
    }
  }else{ 
    cout << "wrong number of arguments" << endl;
  }
}


void TabulatedPotential::Help(string cmd, vector<string> &args)
{
    if(args.size() == 0) {
        if(cmd == "tab") {
            cout << "tab <file> <selection>\n"
                 << "Calculate tabulated potential by inverting the distribution function. "
                    "Statistics is calculated using all interactions in selection.\n"
                    "see also: help tab set\n\n"
                    "example:\ntab set scale bond\ntab U_bond.txt *:bond:*\n";
        }
        if(cmd == "hist") {
            cout << "hist <file> <selection>\n"
                 << "Calculate distribution function for selection. "
                    "Statistics is calculated using all interactions in selection.\n"
                    "see also: help hist set\n\n"
                    "example:hist U_bond.txt *:bond:*\n";
        }
        return;
    }
    if(args[0] == "set") {
        if(args.size() == 1) {
            cout << cmd << " set <option> <value>\n"
                 << "set option for this command. Use \"" << cmd << " set\""
                    " for a list of available options. To get help on a specific option use e.g.\n"
                    << cmd << " set periodic\n";
            return;
        }
        if(args[1] == "n") {
            cout << cmd << "set n <integer>\n"
                 << "set number of bins for table\n";
            return;
        }
        if(args[1] == "min") {
            cout << cmd << "set min <value>\n"
                 << "minimum value of interval for histogram (see also periodic, extend)\n";
            return;
        }
        if(args[1] == "max") {
            cout << cmd << "set max <value>\n"
                 << "maximum value of interval for histogram (see also periodic, extend)\n";
            return;
        }
        if(args[1] == "periodic") {
            cout << cmd << "set periodic <value>\n"
                << "can be 1 for periodic interval (e.g. dihedral) or 0 for "
                    "non-periodic (e.g. bond)\n";
            return;
        }
        if(args[1] == "auto") {
            cout << cmd << "set auto <value>\n"
                 "can be 1 for automatically determine the interval for the "
                 "table (min, max, extend will be ignored) or 0 to use min/max as specified\n";
            return;
        }
        if(args[1] == "extend") {
            cout << cmd << "set extend <value>\n"
                 "should only be used with auto=0. Can be 1 for extend the interval "
                 "if values are out of bounds (min/max) "
                 "or 0 to ignore values which are out of the interal\n";
            return;
        }
        if(args[1] == "scale") {
            cout << cmd << "set scale <value>\n"
                "volume normalization of pdf. Can be no (no scaling), bond "
                "(1/r^2) or angle ( 1/sin(phi) ). See VOTCA manual, section "
                "theoretical background for details\n";
            return;
        }
        if(args[1] == "normalize") {
            cout << cmd << "set normalize <value>\n"
                 "can be 1 for a normalized histogram or 0 to skip normalization\n";
            return;
        }

        if(cmd == "tab") {
            if(args[1] == "smooth_pdf") {
                cout << "tab set smooth_pdf <value>\n"
                 "Perform so many smoothing iterations on the distribution function before inverting the potential\n";
                return;
            }
            if(args[1] == "smooth_pot") {
                cout << "tab set smooth_pot <value>\n"
                 "Perform so many smoothing iterations on tabulated potential after inverting the potential\n";
                return;
            }
            if(args[1] == "T") {
                cout << "tab set T <value>\n"
                 "Temperature in Kelvin the simulation was performed\n";
                return;
            }
        }
    }
    
    cout << "no help text available" << endl;
}

double TabulatedPotential::getTemperature() const { return _Temperature; }

pair<int,int> TabulatedPotential::getSmoothIterations() const {
  return pair<int,int>(_tab_smooth1,_tab_smooth2);
}

void TabulatedPotential::WriteHistogram(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
    DataCollection<double>::selection *sel = nullptr;

    for(size_t i=1; i<args.size(); i++){
        sel = bs.BondedValues().select(args[i], sel);
    }
    Histogram h(_hist_options);
    h.ProcessData(sel);
    out.open(args[0].c_str());
    out << h ;
    out.close();
    cout << "histogram created using " << sel->size() << " data-rows, written to " << args[0] << endl;    
    delete sel;
}

void TabulatedPotential::WritePotential(BondedStatistics &bs, vector<string> &args)
{
    ofstream out;
    DataCollection<double>::selection *sel = nullptr;

    // Appends all the interactions that are specified in args to selection
    // pointer given by sel

    for(size_t i=1; i<args.size(); i++){
      sel = bs.BondedValues().select(args[i], sel);
    }

    Histogram h(_tab_options);        
    
    h.ProcessData(sel);
    for(int i=0; i<_tab_smooth1; ++i){
      Smooth_(h.getPdf(), _tab_options._periodic);
    }
    BoltzmannInvert_(h.getPdf());
    for(int i=0; i<_tab_smooth2; ++i){
        Smooth_(h.getPdf(), _tab_options._periodic);
    }
    out.open(args[0].c_str());
    
    vector<double> F;
    assert(h.getInterval()>0 && "Interval for pdf histogram is 0");
    CalcForce_(h.getPdf(), F, h.getInterval(), _tab_options._periodic);
    for(int i=0; i<h.getN(); i++) {
        out << h.getMin() + h.getInterval()*((double)i) << " " << h.getPdf()[i] << " " << F[i] << endl;
    }
    out.close();
    cout << "histogram created using " << sel->size() << " data-rows, written to " << args[0] << endl;
    delete sel;
}

/******************************************************************************
 * Private Facing Methods
 ******************************************************************************/
bool TabulatedPotential::SetOption_(Histogram::options_t &op, const vector<string> &args)
{
  if(args.size() >2) {
    if(args[1] == "n")
      op._n = lexical_cast<int>(args[2]);       
    else if(args[1] == "min") {
      op._min = lexical_cast<double>(args[2]);
    }else if(args[1] == "max") {
      op._max = lexical_cast<double>(args[2]);
    }else if(args[1] == "periodic"){
      op._periodic = lexical_cast<bool>(args[2]);
    }else if(args[1] == "auto"){
      op._auto_interval = lexical_cast<bool>(args[2]);
    }else if(args[1] == "extend"){
      op._extend_interval = lexical_cast<bool>(args[2]);
    }else if(args[1] == "normalize"){
      op._normalize = lexical_cast<bool>(args[2]);
    }else if(args[1] == "scale") {
      if(args[2]=="no" || args[2]=="bond" || args[2]=="angle"){
        op._scale = args[2];
      }else {
        cout << "scale can be: no, bond or angle\n";
      }
    } else {
      return false;        
    }
  }else {
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

void TabulatedPotential::CalcForce_(vector<double> &U, vector<double> &F, double dx, bool bPeriodic)
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

void TabulatedPotential::Smooth_(vector<double> &data, bool bPeriodic)
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
        old[2] = data[i];
        data[i] = (old[0] + 2.*old[1] + 3.*data[i] + 2.*data[i+1] + data[i+2])/9.;
        old[0]=old[1];
        old[1]=old[2];;
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

void TabulatedPotential::BoltzmannInvert_(vector<double> &data)
{
    double _min, _max;
    
    _min = numeric_limits<double>::max();
    _max = numeric_limits<double>::min();
    
    for(size_t i=0; i<data.size(); i++) {
        _max = max(data[i], _max);
        if(data[i] > 0) _min = min(data[i], _min);
    }
    _max = -kB*_Temperature*log(_max);
    _min = -kB*_Temperature*log(_min)-_max;
    
    for(size_t i=0; i<data.size(); i++) {
        if(data[i] == 0) data[i] = _min;
        else
            data[i] = -kB*_Temperature*log(data[i]) - _max;
    }
}
