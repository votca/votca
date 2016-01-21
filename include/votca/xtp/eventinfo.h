/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/tools/property.h>

namespace votca { namespace xtp {
  
using namespace std;

class Eventinfo {

public:
    
    Eventinfo(){};
    
    void Read_bulk(votca::tools::Property *options){ 

        random_seed                           = options->get("options.general.random_seed").as<int>();
        //const int growsize                    = 2;
        
        electron_density                      = options->get("options.general.electron_density").as<double>();
        hole_density                          = options->get("options.general.hole_density").as<double>();

        number_of_steps                       = options->get("options.general.number_of_steps").as<int>();
        number_of_equilibration_steps         = options->get("options.general.number_of_equilibration_steps").as<int>();

        number_of_report_steps                = options->get("options.general.number_of_report_steps").as<int>();
        
        carrier_trajectory                    = options->get("options.general.carrier_trajectory").as<bool>();
        trajectory_filename                   = options->get("options.general.trajectory_filename").as<string>();
        number_of_trajectory_report_steps     = options->get("options.general.number_of_trajectory_report_steps").as<int>();
        
        filament_visualisation                = options->get("options.general.filament_visualisation").as<bool>();
        visualisation_filename                = options->get("options.general.visualisation_filename").as<string>();
        visualisation_at_nr_steps             = options->get("options.general.visualisation_at_nr_steps").as<int>();

        size_x                                = options->get("options.general.size_x").as<double>();
        size_y                                = options->get("options.general.size_y").as<double>();
        size_z                                = options->get("options.general.size_z").as<double>();
        resize_morphology                     = options->get("options.general.resize_morphology").as<bool>();

        formalism                             = options->get("options.general.formalism").as<string>();        

        alpha                                 = options->get("options.general.alpha").as<double>();
        temperature                           = options->get("options.general.temperature").as<double>();
        efield_x                              = options->get("options.general.efield_x").as<double>();
        efield_y                              = options->get("options.general.efield_y").as<double>();
        efield_z                              = options->get("options.general.efield_z").as<double>();

        gamma                                 = options->get("options.general.gamma").as<double>();
        
        coulomb_strength                      = options->get("options.general.coulomb_strength").as<double>();
        lr_coulomb_strength                   = options->get("options.general.lr_coulomb_strength").as<double>();
        explicit_coulomb                      = options->get("options.general.explicit_coulomb").as<bool>();

        coulomb_cut_off_radius                = options->get("options.general.coulomb_cut_off_radius").as<double>();
        
        binding_energy                        = options->get("options.general.binding_energy").as<double>();
        
        norc                                  = options->get("options.general.norc").as<bool>();
        
        if(advanced) {
            repeat_counting                   = options->get("options.general.repeat_counting").as<bool>();
            no_blocking                       = options->get("options.general.no_blocking").as<bool>();

            mesh_size_x                       = options->get("options.general.mesh_size_x").as<int>();
            mesh_size_y                       = options->get("options.general.mesh_size_y").as<int>();
            mesh_size_z                       = options->get("options.general.mesh_size_z").as<int>();            

            electron_transport_prefactor      = options->get("options.general.electron_prefactor").as<double>();
            hole_transport_prefactor          = options->get("options.general.hole_prefactor").as<double>();
            recombination_prefactor           = options->get("options.general.recombination_prefactor").as<double>();

            visualisation_number_x            = options->get("options.general.visualisation_number_x").as<int>();
            visualisation_number_y            = options->get("options.general.visualisation_number_y").as<int>();
            visualisation_number_z            = options->get("options.general.visualisation_number_z").as<int>();

            number_direct_conv_iv             = options->get("options.general.number_direct_conv_iv").as<int>();
            number_direct_conv_reco           = options->get("options.general.number_direct_conv_reco").as<int>();
        }
        else {
            repeat_counting                   = true;
            no_blocking                       = false;

            mesh_size_x                       = floor(size_x/coulomb_cut_off_radius);
            mesh_size_y                       = floor(size_y/coulomb_cut_off_radius);
            mesh_size_z                       = floor(size_z/coulomb_cut_off_radius); 

            electron_transport_prefactor      = 1.0;
            hole_transport_prefactor          = 1.0;
            recombination_prefactor           = 1.0;

            visualisation_number_x            = floor(size_x);
            visualisation_number_y            = floor(size_y);
            visualisation_number_z            = floor(size_z);
            
            number_direct_conv_iv             = 30;
            number_direct_conv_reco           = 30;
        }
        
    }
    
    void Read_device(votca::tools::Property *options) {
        
        this->Read_bulk(options);
        voltage                               = options->get("options.device.voltage").as<double>();

        init_charges                          = options->get("options.device.init_charges").as<bool>();
        
        left_electron_injection               = options->get("options.device.left_electron_injection").as<bool>();
        left_hole_injection                   = options->get("options.device.left_hole_injection").as<bool>();
        right_electron_injection              = options->get("options.device.right_electron_injection").as<bool>();
        right_hole_injection                  = options->get("options.device.right_hole_injection").as<bool>();

        left_injection_barrier                = options->get("options.device.left_injection_barrier").as<double>();
        right_injection_barrier               = options->get("options.device.right_injection_barrier").as<double>();
        
        left_electrode_distance               = options->get("options.device.left_electrode_distance").as<double>();
        right_electrode_distance              = options->get("options.device.right_electrode_distance").as<double>();

        
        
        timesteps_update_longrange            = options->get("options.device.timesteps_update_longrange").as<int>();
        
        if(advanced) {
            // default value is equal to 1
            injection_prefactor               = options->get("options.device.injection_prefactor").as<double>();
            collection_prefactor              = options->get("options.device.collection_prefactor").as<double>();

            // default value is equal to 0.5
            self_image_prefactor              = options->get("options.device.self_image_prefactor").as<double>();
            
            // default value is integer floor of size of device
            number_of_layers                  = options->get("options.device.number_of_layers").as<int>();
            
            // default value is 10 images
            number_short_range_images         = options->get("options.device.number_short_range_images").as<int>();            
            number_long_range_images          = options->get("options.device.number_long_range_images").as<int>();

            // default slab-based longrange calculations (true)
            longrange_slab                    = options->get("options.device.longrange_slab").as<bool>();

            // default is to write profiles (true))
            write_charge_profile              = options->get("options.device.write_charge_profile").as<bool>();
            write_potential_profile           = options->get("options.device.write_potential_profile").as<bool>();            
           
        }
        else {
            injection_prefactor               = 1.0;
            collection_prefactor              = 1.0;
            
            self_image_prefactor              = 0.5;
            
            number_of_layers                  = floor(size_x);
            
            number_short_range_images         = 10;
            number_long_range_images          = 10;
            
            longrange_slab                    = true;

            write_charge_profile              = true;
            write_potential_profile           = true;      
        }
    }
    
    void Read_cubic(votca::tools::Property *options) {
        NX                                    = options->get("options.cubic.NX").as<int>();
        NY                                    = options->get("options.cubic.NY").as<int>();
        NZ                                    = options->get("options.cubic.NZ").as<int>();
        lat_const                             = options->get("options.cubic.lat_const").as<double>();

        hop_distance                          = options->get("options.cubic.hop_distance").as<double>();
        
        el_disorder                           = options->get("options.cubic.el_disorder").as<double>();
        ho_disorder                           = options->get("options.cubic.ho_disorder").as<double>();
        el_ho_correlation                     = options->get("options.cubic.el_ho_correlation").as<int>();
        
        lumo                                  = options->get("options.cubic.lumo").as<double>();
        homo                                  = options->get("options.cubic.homo").as<double>();

        
    }
    
    void Graph_Parameters(double graph_hopdist, double graph_mindist, votca::tools::vec graph_simboxsize, int graph_maxpairdegree, double av_ho_energy, double stdevhoenergy, double av_elect_energy, double hole_inject_reorg, double electron_inject_reorg) 
    {
        hopdist = graph_hopdist;
        mindist = graph_mindist;
        simboxsize = graph_simboxsize;
        maxpairdegree = graph_maxpairdegree;
        avholeenergy = av_ho_energy;
        stddev_hole_energy = stdevhoenergy;
        avelectronenergy = av_elect_energy;
        hole_injection_reorg = hole_inject_reorg;
        electron_injection_reorg = electron_inject_reorg;

    }
    
    
    void Set_field(){
        efield_z = voltage/(simboxsize.z());
    }
    
    bool advanced; bool device;   
    int random_seed;
    bool init_charges; double hole_density; double electron_density;
    bool explicit_coulomb;

    int number_of_steps; int number_of_equilibration_steps; int number_of_report_steps;
    bool carrier_trajectory; string trajectory_filename; int number_of_trajectory_report_steps;
    bool filament_visualisation; string visualisation_filename; int visualisation_at_nr_steps;
  
    double size_x; double size_y; double size_z; bool resize_morphology;
    string formalism;
    double alpha; double temperature; double efield_x; double efield_y; double efield_z;

    double coulomb_strength; double coulomb_cut_off_radius; double binding_energy; double lr_coulomb_strength;  
    
    int repeat_counting; bool no_blocking;
    int mesh_size_x; int mesh_size_y; int mesh_size_z;
    
    double electron_transport_prefactor; double hole_transport_prefactor; double recombination_prefactor;
    
    int growsize;
    
    int visualisation_number_x; int visualisation_number_y; int visualisation_number_z;
    int number_direct_conv_iv; int number_direct_conv_reco;

    double voltage;
    bool left_electron_injection; bool left_hole_injection; bool right_electron_injection; bool right_hole_injection;
    double left_injection_barrier; double right_injection_barrier; double left_electrode_distance; double right_electrode_distance; 
    int timesteps_update_longrange;
    
    double injection_prefactor; double collection_prefactor;
    double self_image_prefactor;   
    int number_of_layers;
    int number_short_range_images; int number_long_range_images;
    bool longrange_slab;
    bool write_charge_profile; bool write_potential_profile;
           
    double hopdist; double mindist; votca::tools::vec simboxsize; int maxpairdegree; double avholeenergy; double avelectronenergy;
    double stddev_hole_energy;
    double hole_injection_reorg; double electron_injection_reorg;
    
    bool norc;
    
    int NX; int NY; int NZ; double lat_const; double hop_distance;
    double el_disorder; double ho_disorder; int el_ho_correlation;
    double lumo; double homo; double gamma;
};

}} 

#endif

