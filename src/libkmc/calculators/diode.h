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

//#include <votca/kmc/lattice.h>
#include <votca/kmc/state.h>
#include <votca/kmc/graph.h>
#include <votca/kmc/events.h>

using namespace std;

namespace votca { namespace kmc {
    
//typedef votca::tools::vec myvec;

   
class Diode : public KMCCalculator 
{
public:
  Diode() {};
 ~Diode() {};

  void Initialize(Property *options);
  bool EvaluateFrame();

 // State state;


protected:
 void RunKMC(void); 
            
private:
  static const double kB   = 8.617332478E-5; // eV/K
  static const double hbar = 6.5821192815E-16; // eV*s
  static const double eps0 = 8.85418781762E-12/1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
  static const double epsr = 3.0; // relative material permittivity
  static const double Pi   = 3.14159265358979323846;
   
};


void Diode::Initialize(Property *options) {


}

bool Diode::EvaluateFrame()
{
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
}

void Diode::RunKMC() {
    std::cout << "something";
    

}

}}


#endif	/* __VOTCA_KMC_DIODE_H */