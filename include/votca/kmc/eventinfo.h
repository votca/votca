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

#ifndef __VOTCA_KMC_EVENTINFO_H_
#define __VOTCA_KMC_EVENTINFO_H_

namespace votca { namespace kmc {
  
using namespace std;

class Eventinfo {

public:
    
    Eventinfo(){};
    
    void Read(Property *options){ 
    
        nr_electrons                = options->get("options.bulk.nr_electrons").as<int>();
        nr_holes                    = options->get("options.bulk.nr_holes").as<int>();
        nr_reportsteps              = options->get("options.bulk.nr_reportsteps").as<int>();
        traj_store                  = options->get("options.bulk.traj_store").as<bool>();
        traj_filename               = options->get("options.bulk.traj_filename").as<string>();
        
        seed                        = options->get("options.diode.seed").as<int>();
        nr_equilsteps               = options->get("options.diode.nr_equilsteps").as<int>();
        nr_timesteps                = options->get("options.diode.nr_timesteps").as<int>();
        steps_update_longrange      = options->get("options.diode.steps_update_longrange").as<int>();
        
        number_direct_conv_iv       = options->get("options.diode.number_direct_conv_iv").as<int>();
        number_direct_conv_iv       = options->get("options.diode.number_direct_conv_iv").as<int>();
        
        nx                          = options->get("options.diode.nx").as<int>();
        ny                          = options->get("options.diode.ny").as<int>();
        nz                          = options->get("options.diode.nz").as<int>();
        size_x                      = options->get("options.diode.size_x").as<double>();
        size_y                      = options->get("options.diode.size_y").as<double>();
        size_z                      = options->get("options.diode.size_z").as<double>();
        resize                      = options->get("options.diode.resize").as<bool>();
        lattice_constant            = options->get("options.diode.lattice_constant").as<double>();
        left_electrode_distance     = options->get("options.diode.left_electrode_distance").as<double>();
        right_electrode_distance    = options->get("options.diode.right_electrode_distance").as<double>();
        growsize                    = options->get("options.diode.growsize").as<int>();
        alpha                       = options->get("options.diode.alpha").as<double>();
        temperature                 = options->get("options.diode.temperature").as<double>();
        efield_x                    = options->get("options.diode.efield_x").as<double>();
        efield_y                    = options->get("options.diode.efield_y").as<double>();
        efield_z                    = options->get("options.diode.efield_z").as<double>();
        voltage                     = options->get("options.diode.voltage").as<double>();
        injection_barrier           = options->get("options.diode.injection_barrier").as<double>();
        binding_energy              = options->get("options.diode.binding_energy").as<double>();
        formalism                   = options->get("options.diode.formalism").as<string>();
        
        el_density                  = options->get("options.diode.el_density").as<double>();
        ho_density                  = options->get("options.diode.ho_density").as<double>();
        
        electron_prefactor          = options->get("options.diode.electron_prefactor").as<double>();
        hole_prefactor              = options->get("options.diode.hole_prefactor").as<double>();
        injection_prefactor         = options->get("options.diode.injection_prefactor").as<double>();
        collection_prefactor        = options->get("options.diode.collection_prefactor").as<double>();
        recombination_prefactor     = options->get("options.diode.recombination_prefactor").as<double>();
        
        left_electron_injection     = options->get("options.diode.left_electron_injection").as<bool>();
        left_hole_injection         = options->get("options.diode.left_hole_injection").as<bool>();
        right_electron_injection    = options->get("options.diode.right_electron_injection").as<bool>();
        right_hole_injection        = options->get("options.diode.right_hole_injection").as<bool>();
        
        device                      = options->get("options.diode.device").as<int>();
        
        mesh_x                      = options->get("options.diode.mesh_x").as<int>();
        mesh_y                      = options->get("options.diode.mesh_y").as<int>();
        mesh_z                      = options->get("options.diode.mesh_z").as<int>();
        number_layers               = options->get("options.diode.number_layers").as<int>();
        
        nr_sr_images                = options->get("options.diode.nr_sr_images").as<int>();
        nr_lr_images                = options->get("options.diode.nr_lr_images").as<int>();
        
        coulomb_strength            = options->get("options.diode.coulomb_strength").as<double>();
        coulcut                     = options->get("options.diode.coulcut").as<double>();
        self_image_prefactor        = options->get("options.diode.self_image_prefactor").as<double>();
        
        interpolate_longrange       = options->get("options.diode.interpolate_longrange").as<bool>();
        longrange_slab              = options->get("options.diode.longrange_slab").as<bool>();
        
        left_oxide_thickness        = options->get("options.diode.left_oxide_thickness").as<double>();
        right_oxide_thickness       = options->get("options.diode.right_oxide_thickness").as<double>();
        
    }
    
    void Graph_Parameters(double graph_hopdist, double graph_mindist, votca::tools::vec graph_simboxsize, int graph_maxpairdegree, double av_ho_energy, double av_elect_energy) {
        hopdist = graph_hopdist;
        mindist = graph_mindist;
        simboxsize = graph_simboxsize;
        maxpairdegree = graph_maxpairdegree;
        avholeenergy = av_ho_energy;
        avelectronenergy = av_elect_energy;
    }
    
    void Set_field(double in_device_voltage){
        efield_x = voltage/(simboxsize.x());
    }

    int nr_holes; int nr_electrons; int nr_reportsteps; bool traj_store; string traj_filename;
    int seed; int nr_equilsteps; int nr_timesteps; int steps_update_longrange;
    int number_direct_conv_iv; int number_direct_conv_reco;
    int nx; int ny; int nz; int growsize;
    double size_x; double size_y; double size_z; bool resize;
    double lattice_constant; double left_electrode_distance; double right_electrode_distance; 
    double alpha; double temperature; double el_density; double ho_density;
    double efield_x; double efield_y; double efield_z; double voltage;
    double electron_prefactor; double hole_prefactor; double injection_prefactor; double recombination_prefactor; double collection_prefactor;
    string formalism; double injection_barrier; double binding_energy;
    bool left_electron_injection; bool left_hole_injection; bool right_electron_injection; bool right_hole_injection; int device;
    int mesh_x; int mesh_y; int mesh_z; int number_layers;
    double hopdist; double mindist; votca::tools::vec simboxsize; int maxpairdegree; double avholeenergy; double avelectronenergy;
    double coulomb_strength; double coulcut; double self_image_prefactor; int nr_sr_images; long nr_lr_images;
    bool interpolate_longrange; bool longrange_slab;
    double left_oxide_thickness; double right_oxide_thickness;

};

}} 

#endif

