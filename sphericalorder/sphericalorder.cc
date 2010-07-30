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

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include <stdlib.h>
#include <votca/csg/csgapplication.h>
#include <votca/tools/average.h>
#include <votca/tools/tokenizer.h>
#include <votca/csg/cgengine.h>

using namespace std;
using namespace votca::csg;

class CGOrderParam
    : public CsgApplication
{
public:
    
    string ProgramName() {
        return "sphericalorder";
    }

    void HelpText(ostream &out) {

        out << "!! EXPERIMENTAL !! Calculate spherical order parameter.\n"
                " Needs non-spherical beads in mapping.\n\n";

    }

    void Initialize() {
        CsgApplication::Initialize();
        AddProgramOptions()
            ("filter", boost::program_options::value<string>(&_filter)->default_value("*"), "filter molecule names")
            ("radialcut", boost::program_options::value<double>(), "radial cutoff: distance from center where bead is considered")
            ("refmol", boost::program_options::value<string>(&_refmol)->default_value(""), "Reference molecule");
    }

     bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("radialcut");
        return true;
    }

    bool DoTrajectory() { return true;}
    bool DoMapping() { return true;}
    
    void BeginEvaluate(Topology *top, Topology *top_atom) {

        string filter;

        filter = OptionsMap()["filter"].as<string>();
        _radialcutoff = OptionsMap()["radialcut"].as<double>();
        _refmol = OptionsMap()["refmol"].as<string>();

        cout << "using radial cutoff: " << _radialcutoff << endl;

        setFilter(filter);


        _file_u.open("hist_u.xvg");
        if(!_file_u)
                throw runtime_error("cannot open hist_u.xvg for output");
        _file_v.open("hist_v.xvg");
        if(!_file_v)
                throw runtime_error("cannot open hist_v.xvg for output");
        _file_w.open("hist_w.xvg");
        if(!_file_w)
                throw runtime_error("cannot open hist_w.xvg for output");
        
        _n = 0;

        if (_refmol == ""){
            matrix box;
            box = top->getBox();
            vec a = box.getCol(0); vec b = box.getCol(1); vec c = box.getCol(2);
            _ref = (a+b+c)/2;

            cout << "Refernce is center of box " << _ref << endl;
        }
        
        _nbin = 100;
        _hist_u = new real[_nbin];
        _hist_v = new real[_nbin];
        _hist_w = new real[_nbin];

        _nmol=0;

         
    }
    
    void EndEvaluate() {

        cout << "Average number of molecules within cutoff " << (double)_nmol/_n << endl;

        double exp_value = (double)1/_nbin;
        double orderparam = 0;
        
        for (int n=0; n<_nbin; n++){
            _hist_u[n] /= (double)_nmol; // normalize to numberframes and avg. number of molecules
            _hist_v[n] /= (double)_nmol;
            _hist_w[n] /= (double)_nmol;

            _file_u << (double)n*2/(_nbin)-1 << " " << _hist_u[n] << endl;
            _file_v << (double)n*2/(_nbin)-1 << " " << _hist_v[n] << endl;
            _file_w << (double)n*2/(_nbin)-1 << " " << _hist_w[n] << endl;

            orderparam += (_hist_u[n]-exp_value)*(_hist_u[n]-exp_value);
        }

        orderparam = sqrt(orderparam/_nbin);

        cout << "Orderparam " << _radialcutoff << " " << orderparam << endl;

        _file_u.close();
        _file_v.close();
        _file_w.close();
        
    };
    
    void EvalConfiguration(Topology *conf, Topology*conf_atom = 0) {

        vec eR;
        int nmol;
        int nu, nv, nw;
        vec u, v, w;
        nmol=0;

        if (_refmol != ""){
            for(BeadContainer::iterator iter = conf->Beads().begin();
            iter!=conf->Beads().end();++iter) {
                Bead *bead = *iter;
                if(wildcmp("SOL", bead->getName().c_str())){
                    _ref = bead->getPos();
                    cout << " Solute pos " << _ref << endl;
                }

            }
        }
        
        for(BeadContainer::iterator iter = conf->Beads().begin();
        iter!=conf->Beads().end();++iter) {
            Bead *bead = *iter;
            if(!wildcmp(_filter.c_str(), bead->getName().c_str())) continue;
                eR = bead->getPos()-_ref;
                if (abs(eR) < _radialcutoff && abs(eR)>0) {
                    // cout << eR << endl;
                    eR.normalize();
                    u = bead->getU();
                    v = bead->getV();
                    w = bead->getW();
                    u.normalize();
                    v.normalize();
                    w.normalize();

                    nu = (int) ((((eR * u) + 1) / 2) * _nbin);
                    nv = (int) ((((eR * v) + 1) / 2) * _nbin);
                    nw = (int) ((((eR * w) + 1) / 2) * _nbin);

                    //cout << "nu" << nu << "nv" << nv << "nw" << nw << endl;
                    _hist_u[nu] += 1;
                    _hist_v[nv] += 1;
                    _hist_w[nw] += 1;
                    _nmol++;

                }
                
                //cout << "nmol " << _nmol << endl;
        }
        
        _n++;        
    }
    
    void setOut(string filename) { _filename = filename; }
    
    void setFilter(const string &filter) { _filter = filter; }
    
protected:
    ofstream _file;
    string _filename;
    int _n;
    vec  _ref;
    ofstream _file_u;
    ofstream _file_v;
    ofstream _file_w;
    double * _hist_u;
    double * _hist_v;
    double * _hist_w;
    int _nbin;
    int _nmol;
    double _radialcutoff;

    string _filter;
    string _refmol;
};

int main(int argc, char **argv)
{
    CGOrderParam app;
    return app.Exec(argc, argv);
}
