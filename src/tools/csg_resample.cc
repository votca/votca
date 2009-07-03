/* 
 * File:   csg_resample.cc
 * Author: ruehle
 *
 * Created on April 8, 2009, 5:45 PM
 */

#include <tools/cubicspline.h>
#include <tools/table.h>
#include <tools/tokenizer.h>
#include <boost/program_options.hpp>
#include <iostream>

using namespace std;
namespace po = boost::program_options;

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
    Table in, out;

    // program options
    po::options_description desc("Allowed options");            
    
    desc.add_options()
      ("in", po::value<string>(&in_file), "table to read")
      ("out", po::value<string>(&out_file), "table to write")
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
        cout << "csg_resample \n\n";                
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
    for(;out.x(i) < in.x(in.size()-1) && i < out.size(); ++i) {
        while(in.x(j) < out.x(i)) ++j;
        out.flags(i) = in.flags(j);
    }
    
    out.Save(out_file);
    
    return 0;
}

