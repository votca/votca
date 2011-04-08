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

#include <math.h>
#include <votca/tools/histogramnew.h>
#include <votca/tools/tokenizer.h>
#include <csgapplication.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;


class CsgDensityApp
    : public CsgApplication
{
    string ProgramName() { return "csg_density"; }
    void HelpText(ostream &out) { 
        out << "Calculates the mass density distribution along a box axis or radial density profile from reference point";
    }

    // some program options are added here
    void Initialize();
    
    // we want to process a trajectory
    bool DoTrajectory() {return true;}
    bool DoMapping() {return true;}
    bool DoMappingDefault(void) { return false; }
   
    // write out results in EndEvaluate
    void EndEvaluate();
    void BeginEvaluate(Topology *top, Topology *top_atom);
    void EvalConfiguration(Topology *top, Topology *top_ref);

    bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("out", "no output topology specified");
        CheckRequired("trj", "no trajectory file specified");
        return true;
    };

protected:
    string _filter, _out;
    HistogramNew _dist;
    double _rmax;
    int _nbin;
    double _scale;
    int _frames;
    vec _ref;
    vec _axis;
    string _axisname;
    string _molname;
    double _area;
};

int main(int argc, char** argv)
{
    CsgDensityApp app;
    return app.Exec(argc, argv);
}

void CsgDensityApp::BeginEvaluate(Topology *top, Topology *top_atom) {

     matrix box = top->getBox();
     vec a = box.getCol(0);
     vec b = box.getCol(1);
     vec c = box.getCol(2);
    
    _dist.setPeriodic(true);
    _axis=vec(0,0,0);
    _area=0;
    if(_axisname=="x") {
      _axis.setX(1);
      _rmax = abs(a);
      _area= abs(b^c);
    }
    else if(_axisname=="y") {
      _axis.setY(1);
       _rmax = abs(b);
      _area= abs(a^c);
    }
    else if(_axisname=="z") {
      _axis.setZ(1);
      _rmax = abs(c);
      _area= abs(a^b);
    }
    else if(_axisname=="r") {
      _dist.setPeriodic(false);
      _rmax = min(min(abs(a/2), abs(b/2)), abs(c/2));
    } else {
      throw std::runtime_error("unknown axis type");
    }

    if(OptionsMap().count("rmax"))
      _rmax = OptionsMap()["rmax"].as<double>();

    if (_axisname=="r") {
      if(!OptionsMap().count("ref"))
         _ref = a/2+b/2+c/2;
       cout << "Using referece point: " << _ref << endl;
    } 
    else if(OptionsMap().count("ref"))
        throw std::runtime_error("reference center can only be used in case of spherical density");
       
    _dist.Initialize(0, _rmax, _nbin);
    
    cout << "rmax: " << _rmax << endl;
    cout << "axis: " << _axisname << endl;
    cout << "Bins: " << _nbin << endl;
    _frames=0;

}  

void CsgDensityApp::EvalConfiguration(Topology *top, Topology *top_ref)
{
    // loop over all molecules
    bool did_something = false;
    for(MoleculeContainer::iterator imol=top->Molecules().begin(); imol!=top->Molecules().end(); ++imol) {
        Molecule *mol = *imol;
        if(!wildcmp(_molname.c_str(),mol->getName().c_str())) 
            continue;
        int N = mol->BeadCount();
            for(int i=0; i<N; i++) {
                Bead *b = mol->getBead(i);        
                if(!wildcmp(_filter.c_str(), b->getName().c_str()))
                    continue;
                double r;
                if (_axisname=="r") {
                    r = abs(top->BCShortestConnection(_ref, b->getPos()));
                } else {
                    r = b->getPos() *  _axis;
                }
                _dist.Process(r, b->getM());
                did_something = true;
            }
    }
    _frames++;
    if (!did_something) throw std::runtime_error("No molecule in selection");
}


// output everything when processing frames is done
void CsgDensityApp::EndEvaluate()
{
  if (_axisname=="r") {
    _dist.data().y() = _scale/(_frames*_rmax/(double)_nbin *4*M_PI) * element_div( _dist.data().y(),
                   element_prod(_dist.data().x(), _dist.data().x()));
  } else {
    _dist.data().y() = _scale/((double)_frames * _area * _rmax/ (double)_nbin ) *_dist.data().y();
  }
  _dist.data().Save(_out);    
}

// add our user program options
void CsgDensityApp::Initialize()
{
    CsgApplication::Initialize();
    // add program option to pick molecule
    AddProgramOptions("Specific options:")
             ("axis", boost::program_options::value<string>(&_axisname)->default_value("r"), "[x|y|z|r] density axis (r=spherical)")
             ("bins", boost::program_options::value<int>(&_nbin)->default_value(50), "bins")
             ("out", boost::program_options::value<string>(&_out), "Output file")
             ("rmax", boost::program_options::value<double>(), "rmax (default for [r] =min of all box vectors/2, else l )")
             ("scale", boost::program_options::value<double>(&_scale)->default_value(1.0), "scale factor for the density")
             ("molname", boost::program_options::value<string>(&_molname)->default_value("*"), "molname")
             ("filter", boost::program_options::value<string>(&_filter)->default_value("*"), "filter bead names")
             ("ref", boost::program_options::value<vec>(&_ref), "reference zero point");
}

