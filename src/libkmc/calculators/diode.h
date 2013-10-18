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
#include <votca/tools/vec.h>
//#include <votca/kmc/store.h>

using namespace std;

namespace votca { namespace kmc {
    
typedef votca::tools::vec myvec;

enum GraphType{Square, Load};
enum CorrelationType {Correlated, Anticorrelated, Uncorrelated};
enum NodeType {Normal, LeftElectrode, RightElectrode};
    
class Diode : public KMCCalculator 
{
public:
  Diode() {};
 ~Diode() {};

  void Initialize(Property *options);
  bool EvaluateFrame();
  
  int _seed;
  long _equil_time;
  long _run_time;
  long _longrange_update; //number of time steps after which the longrange potential is updated
  
  double _sim_time;
  
  int _injected_holes_left;
  int _injected_holes_right;
  int _injected_electrons_left;
  int _injected_electrons_right;
  int _collected_holes_left;
  int _collected_holes_right;
  int _collected_electrons_left;
  int _collected_electrons_right;
  
  votca::tools::Random2 *RandomVariable;  

  //Lattice _lattice;

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

/*
    if(options->exists("options.lattice.NX")) {
        _lattice->_NX = options->get("options.lattice.NX").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for box size, dimension x. Set to default number (50) \n";
        _lattice->_NX = 50;
    }
    
    if(options->exists("options.lattice.NY")) {
        _lattice->_NY = options->get("options.lattice.NY").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for box size, dimension y. Set to default number (50) \n";
        _lattice->_NY = 50;
    }
    
    if(options->exists("options.lattice.NZ")) {
        _lattice->_NZ = options->get("options.lattice.NZ").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for box size, dimension z. Set to default number (50) \n";
        _lattice->_NZ = 50;
    }
    
    if(options->exists("options.lattice.coul_box_dimX")) {
        _lattice->_coul_box_dimX = options->get("options.lattice.coul_box_dimX").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for coulomb mesh size, dimension x. Set to default number (10) \n";
        _lattice->_coul_box_dimX = 10;
    }

    if(options->exists("options.lattice.coul_box_dimY")) {
        _lattice->_coul_box_dimY = options->get("options.lattice.coul_box_dimY").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for coulomb mesh size, dimension y. Set to default number (10) \n";
        _lattice->_coul_box_dimY = 10;
    }
    
    if(options->exists("options.lattice.coul_box_dimZ")) {
        _lattice->_coul_box_dimZ = options->get("options.lattice.coul_box_dimZ").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for coulomb mesh size, dimension z. Set to default number (10) \n";
        _lattice->_coul_box_dimZ = 10;
    } 

    if(options->exists("options.lattice.site_box_dimX")) {
        _lattice->_site_box_dimX = options->get("options.lattice.site_box_dimX").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for site mesh size, dimension x. Set to default number (2) \n";
        _lattice->_site_box_dimX = 2;
    } 
    
    if(options->exists("options.lattice.site_box_dimY")) {
        _lattice->_site_box_dimY = options->get("options.lattice.site_box_dimY").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for site mesh size, dimension y. Set to default number (2) \n";
        _lattice->_site_box_dimY = 2;
    } 
    
    if(options->exists("options.lattice.site_box_dimZ")) {
        _lattice->_site_box_dimZ = options->get("options.lattice.site_box_dimZ").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for site mesh size, dimension z. Set to default number (2) \n";
        _lattice->_site_box_dimZ = 2;
    } 
    
    if(options->exists("options.lattice.hopping_distance")) {
        _lattice->_hopping_distance= options->get("options.lattice.hopping_distance").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for hopping distance. Set to default number (1) \n";
        _lattice->_hopping_distance = 1;
    } 
    
    if(options->exists("options.lattice.lattice_const")) {
        _lattice->_lattice_const= options->get("options.lattice.lattice_const").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for lattice constant. Set to default number (1) \n";
        _lattice->_lattice_const = 1;
    }    
    
    if(options->exists("options.lattice.alpha")) {
        _lattice->_alpha= options->get("options.lattice.alpha").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for inverse wavefunction decay length. Set to default number (10) \n";
        _lattice->_alpha = 10;
    }
    
    if(options->exists("options.lattice.beta")) {
        _lattice->_beta= options->get("options.lattice.beta").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for 1/(kb T). Set to default number (40) \n";
        _lattice->_beta = 40;
    } 
    
    if(options->exists("options.lattice.efield")) {
        _lattice->_efield= options->get("options.lattice.efield").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for electric field. Set to default number (0.1) \n";
        _lattice->_efield = 0.1;
    }

    if(options->exists("options.lattice.binding_strength")) {
        _lattice->_binding_strength= options->get("options.lattice.binding_strength").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for exciton binding energy. Set to default number (0.3) \n";
        _lattice->_binding_strength = 0.3;
    }     
    
    if(options->exists("options.lattice.coulomb_strength")) {
        _lattice->_coulomb_strength= options->get("options.lattice.coulomb_strength").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for strength of coulomb potential at 1 nm. Set to default number (0.3) \n";
        _lattice->_coulomb_strength = 0.3;
    }     

    if(options->exists("options.lattice.coulomb_cutoff")) {
        _lattice->_coulomb_cutoff= options->get("options.lattice.coulomb_cutoff").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for cut off radius for calculation of coulomb interactions. Set to default number (5) \n";
        _lattice->_coulomb_cutoff = 5;
    }
    
    if(options->exists("options.lattice.electron_prefactor")) {
        _lattice->_electron_prefactor= options->get("options.lattice.electron_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for electron hopping prefactor. Set to default number (1) \n";
        _lattice->_electron_prefactor = 1;
    }
    
    if(options->exists("options.lattice.hole_prefactor")) {
        _lattice->_hole_prefactor= options->get("options.lattice.hole_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for hole hopping prefactor. Set to default number (1) \n";
        _lattice->_hole_prefactor = 1;
    }

    if(options->exists("options.lattice.injection_prefactor")) {
        _lattice->_injection_prefactor= options->get("options.lattice.injection_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for injection prefactor. Set to default number (1) \n";
        _lattice->_injection_prefactor = 1;
    }
 
    if(options->exists("options.lattice.collection_prefactor")) {
        _lattice->_collection_prefactor= options->get("options.lattice.collection_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for collection prefactor. Set to default number (1) \n";
        _lattice->_collection_prefactor = 1;
    }
 
    if(options->exists("options.lattice.recombination_prefactor")) {
        _lattice->_recombination_prefactor= options->get("options.lattice.recombination_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for recombination prefactor. Set to default number (1) \n";
        _lattice->_recombination_prefactor = 1;
    }
    
    if(options->exists("options.lattice.transfer_prefactor")) {
        _lattice->_transfer_prefactor= options->get("options.lattice.transfer_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for normal hopping. Set to default number (1) \n";
        _lattice->_transfer_prefactor = 1;
    }

    if(options->exists("options.lattice.injection_bar_left")) {
        _lattice->_injection_bar_left= options->get("options.lattice.injection_bar_left").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for injection barrier at left electrode. Set to default number (0) \n";
        _lattice->_injection_bar_left = 0;
    }

    if(options->exists("options.lattice.injection_bar_right")) {
        _lattice->_injection_bar_right= options->get("options.lattice.injection_bar_right").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for injection barrier at right electrode. Set to default number (0) \n";
        _lattice->_injection_bar_right = 0;
    }   

    if(options->exists("options.lattice.formalism")) {
        _lattice->_formalism= options->get("options.lattice.formalism").as<string>();
    }
    else {
        cout << "WARNING: Invalid option for hopping formalism. Set to default (Miller Abrahams) \n";
        _lattice->_formalism = Miller;
    }
    
    if(options->exists("options.lattice.correlation")) {
        if(options->get("options.lattice.correlation").as<string>() == "Correlated") {
            _lattice->_correlation = Correlated;
        }
        else if(options->get("options.lattice.correlation").as<string>() == "Anticorrelated") {
            _lattice->_correlation = Anticorrelated;
        }
        else if(options->get("options.lattice.correlation").as<string>() == "Uncorrelated") {
            _lattice->_correlation = Uncorrelated;
        }
        else {
            cout << "Warning: wrong options for energy correlation type (Correlated, Anticorrelated or Uncorrelated). Set to default (Uncorrelated) \n"
            _lattice->_correlation = Uncorrelated;
        }
    }
    else {
        cout << "WARNING: Invalid option for energy correlation. Set to default (Uncorrelated) \n";
        _lattice->_correlation = Uncorrelated;
    }

    if(options->exists("options.lattice.disorder_strength")) {
        _lattice->_disorder_strength= options->get("options.lattice.disorder_strength").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for disorder strength (of holes). Set to default number (0.075) \n";
        _lattice->_disorder_strength = 0.075;
    }

    if(options->exists("options.lattice.disorder_factor")) {
        _lattice->_disorder_factor= options->get("options.lattice.disorder_factor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for ratio of disorder strengths (electron/holes). Set to default number (1.0) \n";
        _lattice->_disorder_factor = 1.0;
    }    

    if(options->exists("options.lattice.graphtype")) {
        if(options->get("options.lattice.graphtype").as<string>() == "Square") {
            _lattice->_graphtype = Square;
        }
        else if(options->get("options.lattice.graphtype").as<string>() == "Load") {
            _lattice->_graphtype = Load;
            
            if(options->exists("options.lattice.graph_filename")) {
              _lattice->_graph_filename= options->get("options.lattice.graph_filename").as<string>();
            }
            else {
              cout << "WARNING: Invalid filename for loading";
            } 

            if(options->exists("options.lattice.load_pairs")) {
              _lattice->_load_pairs= options->get("options.lattice.load_pairs").as<bool>();
            }
            else {
              cout << "WARNING: Invalid option for loading of pairs or not. Set to default number (True) \n";
              _lattice->_load_pairs = true;
            } 

            if(options->exists("options.lattice.load_static_energies")) {
              _lattice->_load_static_energies= options->get("options.lattice.load_static_energies").as<bool>();
            }
            else {
              cout << "WARNING: Invalid option for loading of static energies or not. Set to default number (True) \n";
              _lattice->_load_static_energies = true;
            } 

            if(options->exists("options.lattice.load_transfers")) {
              _lattice->_load_transfers= options->get("options.lattice.load_transfers").as<bool>();
            }
            else {
              cout << "WARNING: Invalid option for loading of transfer integrals or not. Set to default number (True) \n";
              _lattice->_load_transfers = true;
            }
            
        }
        else {
            cout << "Warning: wrong options for graph type (Square or Load). Set to default (Square) \n"
            _lattice->_graphtype = Square;
        }
    }
    else {
        cout << "WARNING: Invalid option for graph type. Set to default (Square) \n";
        _lattice->_graphtype = Square;
    }
    
    if(options->exists("options.lattice.init_maxholes")) {
        _lattice->_init_maxholes= options->get("options.lattice.init_maxholes").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for initial number of maximum supported number of holes. Set to default number (0) \n";
        _lattice->_init_maxholes= 0;
    } 

    if(options->exists("options.lattice.init_maxelectrons")) {
        _lattice->_init_maxelectrons= options->get("options.lattice.init_maxelectrons").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for initial number of maximum supported number of electrons. Set to default number (0) \n";
        _lattice->_init_maxelectrons= 0;
    }

    if(options->exists("options.lattice.init_maxexcitons")) {
        _lattice->_init_maxexcitons= options->get("options.lattice.init_maxexcitons").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for initial number of maximum supported number of excitons. Set to default number (0) \n";
        _lattice->_init_maxexcitons= 0;
    }

    if(options->exists("options.lattice.number_of_shortrange_images")) {
        _lattice->_number_of_shortrange_images= options->get("options.lattice.number_of_shortrange_images").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for number of images in short range coulomb potential calculations. Set to default number (1) \n";
        _lattice->_number_of_shortrange_images = 1;
    } 

    if(options->exists("options.lattice.number_of_longrange_images")) {
        _lattice->_number_of_longrange_images= options->get("options.lattice.number_of_longrange_images").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for number of images in long range coulomb potential calculations. Set to default number (100000) \n";
        _lattice->_number_of_longrange_images = 100000;
    }     
    
    if(options->exists("options.lattice.left_electrode_distance")) {
        _lattice->_left_electrode_distance= options->get("options.lattice.left_electrode_distance").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for distance of lattice to left electrode. Set to default number (1.0) \n";
        _lattice->_left_electrode_distance = 1.0;
    }

    if(options->exists("options.lattice.right_electrode_distance")) {
        _lattice->_right_electrode_distance= options->get("options.lattice.right_electrode_distance").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for distance of lattice to right electrode. Set to default number (1.0) \n";
        _lattice->_right_electrode_distance = 1.0;
    }   

    if(options->exists("options.lattice.device")) {
        _lattice->_device= options->get("options.lattice.device").as<bool>();
    }
    else {
        cout << "WARNING: Invalid option for device setup or not?. Set to default value of true \n";
        _lattice->_device = true;
    }
    
    if(options->exists("options.lattice.dual_injection")) {
        _lattice->_dual_injection= options->get("options.lattice.dual_injection").as<bool>();
    }
    else {
        cout << "WARNING: Invalid option for injection of both holes and electrons at both electrodes at the same time?. Set to default value of false \n";
        _lattice->_dual_injection = false;
    }

    if(options->exists("options.lattice.self_image_prefactor")) {
        _lattice->_self_image_prefactor = options->get("options.lattice.self_image_prefactor").as<double>();
    }
    else {
        cout << "WARNING: Invalid option for self image prefactor. Set to default number (0.5) \n";
        _lattice->_self_image_prefactor = 0.5;
    }
    
    if(options->exists("options.lattice.seed")) {
        _seed = options->get("options.lattice.seed").as<int>();
    }
    else {
        cout << "WARNING: Invalid option for random seed. Set to default number (1) \n";
        _seed = 1;
    }
    
    if(options->exists("options.lattice.equil_time")) {
        _equil_time = options->get("options.lattice.equil_time").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for number of equilibration steps. Set to default number (50000) \n";
        _equil_time = 50000;
    }
    
    if(options->exists("options.lattice.run_time")) {
        _run_time = options->get("options.lattice.run_time").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for number of simulation steps. Set to default number (100000) \n";
        _run_time = 100000;
    }
    
    if(options->exists("options.lattice.longrange_update")) {
        _longrange_update = options->get("options.lattice.longrange_update").as<long>();
    }
    else {
        cout << "WARNING: Invalid option for number after which the longrange coulomb potential is updated. Set to default number (100) \n";
        _longrange_update = 100;
    }

    _
    srand(_seed); // srand expects any integer in order to initialise the random number generator
    RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());

    _lattice->Initialize(RandomVariable);
    _lattice->_longrange->Reset();
    
    _sim_time = 0.0;
 */
}

bool Diode::EvaluateFrame()
{
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
}

void Diode::RunKMC() {

    /*
    // Loop over timesteps:
    for (long itstep=0; itstep<2*NEQUIL+NTIMESTEPS; itstep++) {
      if (div(itstep,_longrange_update).rem==0 && itstep>0) {
        // Update longrange-interactions
        _lattice->_longrange->Update_cache();
        _lattice->Recompute_Injection_Rates();
        _lattice->Recompute_Jump_Rates();
      }
      
      //obtain random number for vssm group

      double electrontotalrate = _lattice->_electron_rates.recompute();
      double holetotalrate = _lattice->_hole_rates.recompute();
      double totalrate = electrontotalrate + holetotalrate;
      
      double randomeorh = RandomVariable->rand_uniform(totalrate);

       // Find chosen charge, chosen jump:
      int chosencharge;
      int chosenjump;
      int charge;
      Carriertype cartype1;

      Node* initnode;
      Node* finalnode;      
      
      if (randomeorh < electrontotalrate) { //electron
        double randrate = randomeorh;
        _lattice->_electron_rates.search(chosencharge,chosenjump,randrate); // Note: chosencharge and chosenjump are passed by reference
        charge = -1;
        cartype1 = Electron;
        if(chosencharge == -1) {
            initnode = _lattice->_graph->_right_electrode;
        }
        else {
            int node_index = _lattice->_electron_state.Get_item(chosencharge)->_node_index;
            initnode = _lattice->_graph->_nodes[_node_index];
        }
      }
      else {
        double randrate = randomeorh - electrontotalrate;
        _lattice->_hole_rates.search(chosencharge,chosenjump,randrate);
        charge = 1;
        cartype1 = Hole;
        if(chosencharge == -1) {
            initnode = _lattice->_graph->_left_electrode;
        }
        else {
            int node_index = _lattice->_hole_state.Get_item(chosencharge)->_node_index;
            initnode = _lattice->_graph->_nodes[_node_index];
        }
      }

      // Update simulation time
      double rate_inv = 1.0/totalrate;
      double randomtime = 1- RandomVariable->rand_uniform(totalrate);
      double deltatime = (-1.0*rate_inv)*log(randomtime);    
      _sim_time += deltatime;


      // Execute the chosen event and update rates
      if (chosencharge==-1) { // Injection event (and possible recombination)
          if(charge==1) {
              finalnode = initnode->_ho_pair_nodes[chosenjump];
              _injected_electrons_right++;
//              LayerCurrent.adjust(NX,1);
          }
          else {
              finalnode = initnode->_el_pair_nodes[chosenjump];
              _injected_holes_left++;
//              LayerCurrent.adjust(0,1);
          }
          _downfield_hops++;
          int my_elcarrier = _lattice->finalnode->_electron_number;
          int my_hocarrier = _lattice->finalnode->_hole_number;
          if ((my_elcarrier == -1)&&(my_hocarrier == -1)) { // Injection to empty site
              _lattice->Add_Remove_Carrier(finalnode,cartype1,Add);
          }
          else if(my_elcarrier != -1) {
              if(cartype1 == Electron) { //Blocking
              }
              else if(cartype1 == Hole) { // Recombination
                _lattice->Add_Remove_Carrier(finalnode,Electron,Remove);
//              recombinations++;
//              RecombinationDensity.adjust(s2.x,1);
              }
          }
          else if(my_hocarrier != -1) {
              if(cartype1 == Hole) { //Blocking (need to recognize bipolarons?)
              }
              else if(cartype1 == Electron) {
                  _lattice->Add_Remove_Carrier(finalnode,Hole,Remove);
//              recombinations++;
//              RecombinationDensity.adjust(s2.x,1);
              }
          }
      }
      else {
          if(charge==1) {
              finalnode  = initnode->_ho_pair_nodes[chosenjump];
          }
          else if(charge==-1) {
              finalnode = initnode->_el_pair_nodes[chosenjump];
          }
//        LayerCurrent.adjust(s1.x+(my_jump.x+1)/2,charge*my_jump.x);
          if (finalnode->_node_type != Normal) { // Collection event
            if (finalnode->_node_type == RightElectrode) { // Collection event at cathode
              if (charge==1) {
                _collected_holes_right++;
              }
              else {
                _collected_electrons_right++;
              }
//              downfield_hops += charge;
            }
            else { // Collection event at anode
              if (charge==1) {
                _collected_holes_left++;
              }
              else {
                _collected_electrons_left++;
              }
//              downfield_hops -= charge;
            }
            my_lattice->add_remove_carrier(s1,cartype1,Remove);

          }
          else { // Not a collection event
//          downfield_hops += charge*my_jump.x;
            
            int my_elcarrier = _lattice->finalnode->_electron_number;
            int my_hocarrier = _lattice->finalnode->_hole_number;             
              
            if (((my_elcarrier != -1)&&(cartype1==Hole))||((my_hocarrier != -1)&&(cartype1==Electron))) { // Recombination event
              my_lattice->add_remove_carrier(initnode,cartype1,Remove);
              Carriertype cartype2
              if(cartype1 == Hole) {
                  cartype2 = Electron;
              }
              else {
                  cartype2 = Hole;
              }
              my_lattice->add_remove_carrier(finalnode,cartype2,Remove);
//              recombinations++;
//              RecombinationDensity.adjust(s2.x,1);
            }
            else { // Ordinary hop (blocking should be included)
              my_lattice->add_remove_carrier(initnode,cartype1,Remove);
              my_lattice->add_remove_carrier(finalnode,cartype,Add);
            }
          }
        }
    }














    
*/

}

}}


#endif	/* __VOTCA_KMC_DIODE_H */