/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_DIODE_H
#define	__VOTCA_KMC_DIODE_H

#include <votca/kmc/graph.h>
#include <votca/tools/vec.h>
#include <votca/kmc/carrier.h>
#include <votca/kmc/state.h>
#include <votca/kmc/event.h>

using namespace std;

namespace votca { namespace kmc {
    
class Diode : public KMCCalculator 
{
public:
  Diode() {};
 ~Diode() {};

  // void Initialize(const char *filename, Property *options, const char *outputfile );
  void Initialize(Property *options);
  bool EvaluateFrame();
  
  void RunKMC(void);
            
protected:
    
 string _lattice_type;
 
 // square lattice
 
 int _Nbox_x;
 int _Nbox_y;
 int _Nbox_z;
 double _lattice_const;
            
private:
  static const double kB   = 8.617332478E-5; // eV/K
  static const double hbar = 6.5821192815E-16; // eV*s
  static const double eps0 = 8.85418781762E-12/1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
  static const double epsr = 3.0; // relative material permittivity
  static const double Pi   = 3.14159265358979323846;
  
  //input/output/options file
  //string _statefile;
  //string _optionsxml;
  //string _outputfile;
  
    Graph _graph;
   
};


// void Diode::Initialize(const char *filename, Property *options, const char *outputfile ) {
void Diode::Initialize(Property *options) {   
        if (options->exists("options.lattice.lattice_type")) {
	    _lattice_type = options->get("options.lattice.lattice_type").as<string>();
	}
        else {
	    cout << "Reading node coordinates from state file \n";
            _lattice_type = "statefile";
        }
        if(_lattice_type != "square" && _lattice_type != "statefile"){
	    cout << "WARNING in simulation: Invalid options for lattice type. Correct options: statefile, square. Set to default option (statefile) \n" ;
            _lattice_type = "statefile";
        }

        if(_lattice_type == "square") {
            if(options->exists("options.lattice.Nbox_x")) {
                 _Nbox_x = options->get("options.lattice.Nbox_x").as<int>();
            }
            else {
                 cout << "WARNING: Invalid option for box size, dimension x. Set to default number (50 sites) \n";
                 _Nbox_x = 50;
            }
            if(options->exists("options.lattice.Nbox_y")) {
                 _Nbox_y = options->get("options.lattice.Nbox_y").as<int>();
            }
            else {
                 cout << "WARNING: Invalid option for box size, dimension y. Set to default number (50 sites) \n";
                 _Nbox_y = 50;
            }
            if(options->exists("options.lattice.Nbox_z")) {
                 _Nbox_z = options->get("options.lattice.Nbox_z").as<int>();
            }
            else {
                 cout << "WARNING: Invalid option for box size, dimension z. Set to default number (50 sites) \n";
                 _Nbox_z = 50;
            }            
            if(options->exists("options.lattice.lattice_const")) {
                 _lattice_const = options->get("options.lattice.lattice_const").as<double>();
            }
            else {
                 cout << "WARNING: Invalid option for lattice constant. Set to default number (1.0) \n";
                 _lattice_const = 1.0;
            }
        }
        

//        _filename = filename;
//        _outputfile = outputfile;    
   
    
}

bool Diode::EvaluateFrame()
{
    RunKMC();
}

void Diode::RunKMC() {
    
    cout << "I am in Run KMC\n" ;
 //   cout << _graph->nodes[0]->nodeposition.x;


    _graph.CreateSquareLattice(_Nbox_x,_Nbox_y,_Nbox_z,_lattice_const);
    
    State _state;
    _state.Load();

}

}}


#endif	/* __VOTCA_KMC_DIODE_H */