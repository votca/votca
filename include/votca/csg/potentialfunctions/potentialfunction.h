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

#ifndef POTENTIALFUNCTION_H
#define	POTENTIALFUNCTION_H

#include <votca/tools/table.h>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>

using namespace std;
using namespace votca::tools;

class PotentialFunction {
public:

    virtual ~PotentialFunction() {}
    // read parameters from the input file
    virtual void setParam(string filename);
    // save parameters to the file
    virtual void SaveParam(const string& filename);
    // write potential table
    virtual void SavePotTab(const string& filename, const double step);
    // write potential table for specified interval
    virtual void SavePotTab(const string& filename, const double step,
			    const double rmin, const double rcut);
    // set all parameters
    void setParam(const Eigen::VectorXd& param){ _lam = param; }
    // set ith parameter
    void setParam(const int i, const double val) { _lam(i) = val; }
    // set ith parameter among those to be optimized
    virtual void setOptParam(const int i, const double val) {
        setParam(i,val);
    }
    // set minimum r value to avoid large values
    void setMinDist(const double min) { _min = min; }
    // set cut-off value
    void setCutOffDist(const double cutoff) { _cut_off = cutoff; }
    // calculate function
    virtual double CalculateF (const double r) const = 0;
    // calculate first derivative w.r.t. ith parameter
    virtual double CalculateDF(const int i, const double r) const = 0;
    // calculate second derivative w.r.t. ith parameter
    virtual double CalculateD2F(const int i, const int j, const double r) const = 0;
    // return parameter
    Eigen::VectorXd& Params() { return _lam; }
    // return ith parameter
    double getParam(const int i) const { return _lam(i); }
    // return ith parameter among those to be optimized
    virtual double getOptParam(const int i) const {
        return getParam(i);
    }
    // return size of parameters
    int getParamSize() const { return _lam.size(); }
    // return size of parameters to be optimized
    virtual int getOptParamSize() const { return getParamSize();}
    // return cut-off value
    double getCutOff() const { return _cut_off; }
    double getMinDist() const { return _min; }

protected:

    PotentialFunction(const string& name_,const int nlam_,const double min_,const double max_);

    string _name;
    Eigen::VectorXd _lam;
    double _cut_off;
    double _min;
};

#endif
