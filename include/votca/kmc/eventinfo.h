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
    
    void Read_bulk(Property *options){ 

        seed                        = options->get("options.general.seed").as<int>();
        
        init_charges                = options->get("options.general.init_charges").as<bool>();
        int_charge_readout          = options->get("options.general.int_charge_readout").as<bool>();
        nr_electrons                = options->get("options.general.nr_electrons").as<int>();
        nr_holes                    = options->get("options.general.nr_holes").as<int>();
        el_density                  = options->get("options.general.el_density").as<double>();
        ho_density                  = options->get("options.general.ho_density").as<double>();

        nr_timesteps                = options->get("options.general.nr_timesteps").as<int>();
        nr_equilsteps               = options->get("options.general.nr_equilsteps").as<int>();
        nr_reportsteps              = options->get("options.general.nr_reportsteps").as<int>();
        repeat_counting             = options->get("options.general.repeat_counting").as<bool>();
        
        number_direct_conv_iv       = options->get("options.general.number_direct_conv_iv").as<int>();
        number_direct_conv_reco     = options->get("options.general.number_direct_conv_reco").as<int>();

        nr_traj_reportsteps         = options->get("options.general.nr_traj_reportsteps").as<int>();
        traj_store                  = options->get("options.general.traj_store").as<bool>();
        traj_filename               = options->get("options.general.traj_filename").as<string>();
        
        viz_store                   = options->get("options.general.viz_store").as<bool>();
        viz_filename                = options->get("options.general.viz_filename").as<string>();
        viz_nr_timesteps            = options->get("options.general.viz_nr_timesteps").as<int>();
        viz_nx                      = options->get("options.general.viz_nx").as<int>();
        viz_ny                      = options->get("options.general.viz_ny").as<int>();
        viz_nz                      = options->get("options.general.viz_nz").as<int>();

        size_x                      = options->get("options.general.size_x").as<double>();
        size_y                      = options->get("options.general.size_y").as<double>();
        size_z                      = options->get("options.general.size_z").as<double>();
        resize                      = options->get("options.general.resize").as<bool>();

        rate_calculate              = options->get("options.general.rate_calculate").as<bool>();
        formalism                   = options->get("options.general.formalism").as<string>();        
        growsize                    = options->get("options.general.growsize").as<int>();

        alpha                       = options->get("options.general.alpha").as<double>();
        temperature                 = options->get("options.general.temperature").as<double>();
        efield_x                    = options->get("options.general.efield_x").as<double>();
        efield_y                    = options->get("options.general.efield_y").as<double>();
        efield_z                    = options->get("options.general.efield_z").as<double>();
        binding_energy              = options->get("options.general.binding_energy").as<double>();
        electron_prefactor          = options->get("options.general.electron_prefactor").as<double>();
        hole_prefactor              = options->get("options.general.hole_prefactor").as<double>();
        recombination_prefactor     = options->get("options.general.recombination_prefactor").as<double>();

        mesh_x                      = options->get("options.general.mesh_x").as<int>();
        mesh_y                      = options->get("options.general.mesh_y").as<int>();
        mesh_z                      = options->get("options.general.mesh_z").as<int>();
        coulomb_strength            = options->get("options.general.coulomb_strength").as<double>();
        coulcut                     = options->get("options.general.coulcut").as<double>();
        
        no_blocking                 = options->get("options.general.no_blocking").as<bool>();
        device                      = options->get("options.device.device").as<bool>();

        fixed_car                   = options->get("options.general.fixed_car").as<bool>();        

    }
    
    void Read_device(Property *options) {
        
        this->Read_bulk(options);
        voltage                     = options->get("options.device.voltage").as<double>();
        injection_prefactor         = options->get("options.device.injection_prefactor").as<double>();
        collection_prefactor        = options->get("options.device.collection_prefactor").as<double>();
        left_electron_injection     = options->get("options.device.left_electron_injection").as<bool>();
        left_hole_injection         = options->get("options.device.left_hole_injection").as<bool>();
        right_electron_injection    = options->get("options.device.right_electron_injection").as<bool>();
        right_hole_injection        = options->get("options.device.right_hole_injection").as<bool>();

        left_injection_barrier      = options->get("options.device.left_injection_barrier").as<double>();
        right_injection_barrier     = options->get("options.device.right_injection_barrier").as<double>();
        left_electrode_distance     = options->get("options.device.left_electrode_distance").as<double>();
        right_electrode_distance    = options->get("options.device.right_electrode_distance").as<double>();
        
        hole_inject_reorg           = options->get("options.device.hole_inject_reorg").as<double>();
        electron_inject_reorg        = options->get("options.device.electron_inject_reorg").as<double>();
        
        number_layers               = options->get("options.device.number_layers").as<int>();
        nr_sr_images                = options->get("options.device.nr_sr_images").as<int>();
        self_image_prefactor        = options->get("options.device.self_image_prefactor").as<double>();
        lr_adaptive                 = options->get("options.device.lr_adaptive").as<bool>();
        nr_lr_images                = options->get("options.device.nr_lr_images").as<int>();
        lr_images_accuracy          = options->get("options.device.lr_images_accuracy").as<double>();        
        interpolate_longrange       = options->get("options.device.interpolate_longrange").as<bool>();
        longrange_slab              = options->get("options.device.longrange_slab").as<bool>();
        steps_update_longrange      = options->get("options.device.steps_update_longrange").as<int>();
        write_charge_profile        = options->get("options.device.write_charge_profile").as<bool>();
        write_pot_profile           = options->get("options.device.write_pot_profile").as<bool>();
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
    int viz_nx; int viz_ny; int viz_nz; int growsize;
    double size_x; double size_y; double size_z; bool resize;
    double lattice_constant; double left_electrode_distance; double right_electrode_distance; 
    double alpha; double temperature; double el_density; double ho_density;
    double efield_x; double efield_y; double efield_z; double voltage;
    double electron_prefactor; double hole_prefactor; double injection_prefactor; double recombination_prefactor; double collection_prefactor;
    string formalism; double left_injection_barrier; double right_injection_barrier; double binding_energy;
    bool left_electron_injection; bool left_hole_injection; bool right_electron_injection; bool right_hole_injection; int device;
    int mesh_x; int mesh_y; int mesh_z; int number_layers;
    double hopdist; double mindist; votca::tools::vec simboxsize; int maxpairdegree; double avholeenergy; double avelectronenergy;
    double coulomb_strength; double coulcut; double self_image_prefactor; int nr_sr_images; long nr_lr_images;
    bool interpolate_longrange; bool longrange_slab;
    bool rate_calculate; bool viz_store; string viz_filename; int viz_nr_timesteps;
    bool lr_adaptive; double lr_images_accuracy;
    bool int_charge_readout; int nr_traj_reportsteps; bool repeat_counting;
    double hole_inject_reorg; double electron_inject_reorg;
    bool no_blocking; bool init_charges; bool write_charge_profile; bool write_pot_profile; bool fixed_car;

};

}} 

#endif

