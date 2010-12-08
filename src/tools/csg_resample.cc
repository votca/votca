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
#include <votca/tools/akimaspline.h>
#include <votca/tools/linspline.h>
#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>
#include <boost/program_options.hpp>
#include <iostream>
#include "version.h"

using namespace std;
namespace po = boost::program_options;
using namespace votca::csg;
using namespace votca::tools;

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

    string in_file, out_file, grid, fitgrid, comment, type, boundaries;
    CubicSpline spline;
    AkimaSpline akspline;
    LinSpline linspline;
    Table in, out, der;

    // program options
    po::options_description desc("Allowed options");            
    
    desc.add_options()
      ("in", po::value<string>(&in_file), "table to read")
      ("out", po::value<string>(&out_file), "table to write")
      ("derivative", po::value<string>(), "table to write")
      ("grid", po::value<string>(&grid), "new grid spacing (min:step:max). If 'grid' is specified only, interpolation is performed.")
      ("type", po::value<string>(&type)->default_value("akima"), "[cubic|akima|linear]. If option is not specified, the default type 'akima' is assumed.")
      ("fitgrid", po::value<string>(&fitgrid), "specify fit grid (min:step:max). If 'grid' and 'fitgrid' are specified, a fit is performed.")
      ("comment", po::value<string>(&comment), "store a comment in the output table")
      ("boundaries", po::value<string>(&boundaries), "(natural|periodic|derivativezero) sets boundary conditions")
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

    if(!(vm.count("grid") || vm.count("fitgrid"))) {
            cout << "Need grid for interpolation or fitgrid for fit.\n";
            return 1;
    }

    if((!vm.count("grid")) && vm.count("fitgrid")) {
            cout << "Need a grid for fitting as well.\n";
            return 1;
    }

    
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

    if (vm.count("type")) {
        if(type=="cubic") {
            spline.setBC(CubicSpline::splineNormal);
        }
        if(type=="akima") {
            akspline.setBC(AkimaSpline::splineNormal);
        }
        if(type=="linear") {
            linspline.setBC(LinSpline::splineNormal);
        }
    }
    

    if (vm.count("boundaries")) {
        if(boundaries=="periodic") {
            spline.setBC(CubicSpline::splinePeriodic);
            akspline.setBC(AkimaSpline::splinePeriodic);
            linspline.setBC(LinSpline::splinePeriodic);
        }
        if(boundaries=="derivativezero") {
            spline.setBC(CubicSpline::splineDerivativeZero);
        }
        //default: normal
    }    


    // in case fit is specified
    if (vm.count("fitgrid")) {
        Tokenizer tok(fitgrid, ":");
        vector<string> toks;
        tok.ToVector(toks);
        if(toks.size()!=3) {
            cout << "wrong range format in fitgrid, use min:step:max\n";
            return 1;        
        }
        double sp_min, sp_max, sp_step;
        sp_min = boost::lexical_cast<double>(toks[0]);
        sp_step = boost::lexical_cast<double>(toks[1]);
        sp_max = boost::lexical_cast<double>(toks[2]);
        cout << "doing " << type << " fit " << sp_min << ":" << sp_step << ":" << sp_max << endl;

        if(type=="cubic") {
            spline.GenerateGrid(sp_min, sp_max, sp_step);
            spline.Fit(in.x(), in.y());
        }
        if(type=="akima") {
            akspline.GenerateGrid(sp_min, sp_max, sp_step);
            akspline.Fit(in.x(), in.y());
        }
        if(type=="linear") {
            linspline.GenerateGrid(sp_min, sp_max, sp_step);
            linspline.Fit(in.x(), in.y());
        }
    } else {
        // else: do interpolation (default = cubic)
        if(type=="cubic") {
            spline.Interpolate(in.x(), in.y());
        }
        if(type=="akima") {
            akspline.Interpolate(in.x(), in.y());
        }
        if(type=="linear") {
            linspline.Interpolate(in.x(), in.y());
        }
    }

    
    out.GenerateGridSpacing(min, max, step);

    if(type=="cubic") {
        spline.Calculate(out.x(), out.y());
    }
    if(type=="akima") {
        akspline.Calculate(out.x(), out.y());
    }
    if(type=="linear") {
        linspline.Calculate(out.x(), out.y());
    }
    
    
    //store a comment line
    if (vm.count("comment")){
        out.set_comment(comment);
    }
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

        if (type == "cubic") {
            spline.CalculateDerivative(der.x(), der.y());
        }
        if (type == "akima") {
            akspline.CalculateDerivative(der.x(), der.y());
        }
        if (type == "linear") {
            linspline.CalculateDerivative(der.x(), der.y());
        }
        der.Save(vm["derivative"].as<string>());
    }
    return 0;
}

