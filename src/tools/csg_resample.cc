/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/tools/cubicspline.h>
#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>
#include <boost/program_options.hpp>
#include <iostream>
#include "version.h"

using namespace std;
namespace po = boost::program_options;

void help_text()
{
    votca::csg::HelpTextHeader("csg_resample");
    cout << "Change grid + interval of any sort of table files.\n"
            "Mainly called internally by inverse script, can also be\n"
            "used to manually prepare input files for coarse-grained\n"
            "simulations.\n\n";     
}

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
    if(!vm.count(option)) {
        cout << "csg_resample \n\n";                
        cout << desc << endl << "parameter " << option << " is not specified\n";
        exit(1);
    }
}

int main(int argc, char** argv)
{
    string in_file, out_file, grid, spfit;
    CubicSpline spline;
    Table in, out, der;

    // program options
    po::options_description desc("Allowed options");            
    
    desc.add_options()
      ("in", po::value<string>(&in_file), "table to read")
      ("out", po::value<string>(&out_file), "table to write")
      ("derivative", po::value<string>(), "table to write")
      ("grid", po::value<string>(&grid), "new grid spacing (min:step:max)")
      ("spfit", po::value<string>(&spfit), "specify spline fit grid. if option is not specified, normal spline interpolation is performed")
      //("bc", po::)
      ("help", "options file for coarse graining");
    
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }
    
    // does the user want help?
    if (vm.count("help")) {
        help_text();
        cout << desc << endl;
        return 0;
    }
    
    check_option(desc, vm, "in");
    check_option(desc, vm, "out");
    check_option(desc, vm, "grid");
    
    double min, max, step;
    {
        Tokenizer tok(grid, ":");
        vector<string> toks;
        tok.ToVector(toks);
        if(toks.size()!=3) {
            cout << "wrong range format, use min:step:max\n";
            return 1;        
        }
        min = boost::lexical_cast<double>(toks[0]);
        step = boost::lexical_cast<double>(toks[1]);
        max = boost::lexical_cast<double>(toks[2]);
    }
        
    
    in.Load(in_file);
    spline.setBC(CubicSpline::splineNormal);

    if (vm.count("spfit")) {
        Tokenizer tok(spfit, ":");
        vector<string> toks;
        tok.ToVector(toks);
        if(toks.size()!=3) {
            cout << "wrong range format in spfit, use min:step:max\n";
            return 1;        
        }
        double sp_min, sp_max, sp_step;
        sp_min = boost::lexical_cast<double>(toks[0]);
        sp_step = boost::lexical_cast<double>(toks[1]);
        sp_max = boost::lexical_cast<double>(toks[2]);
        cout << "doing spline fit " << sp_min << ":" << sp_step << ":" << sp_max << endl;
        spline.GenerateGrid(sp_min, sp_max, sp_step);
                
        spline.Fit(in.x(), in.y());
    } else {
        spline.Interpolate(in.x(), in.y());
    }
    
    out.GenerateGridSpacing(min, max, step);
    spline.Calculate(out.x(), out.y());
    
    out.y() = out.y();
    out.flags() = ub::scalar_vector<double>(out.flags().size(), 'o');

    int i=0;
    for(i=0; out.x(i) < in.x(0) && i<out.size(); ++i);

    int j=0;
    for(;i < out.size(); ++i) {
        for(; j < in.size(); ++j)
            if(in.x(j) >= out.x(i))
                break;        
        if(in.size() == j) break;
        out.flags(i) = in.flags(j);
    }
    
    out.Save(out_file);
    
    if (vm.count("derivative")) {
        der.GenerateGridSpacing(min, max, step);
        der.flags() = ub::scalar_vector<double>(der.flags().size(), 'o');
        spline.CalculateDerivative(der.x(), der.y());
        der.Save(vm["derivative"].as<string>());
    }
    return 0;
}

