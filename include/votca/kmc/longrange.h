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

#ifndef __VOTCA_KMC_LONGRANGE_H_
#define __VOTCA_KMC_LONGRANGE_H_

#include <votca/tools/vec.h>
#include <votca/kmc/graph.h>
#include <votca/kmc/globaleventinfo.h>

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;

class Longrange {

public:
    void Update_cache(myvec sim_box_size, Globaleventinfo* globevent); // Update cached longrange contributions
    double Get_cached_longrange(int layer); // Return cached value
    void Reset();
    //note that the number of images for the calculation of the long range potential should be considerably larger 
    //than the number for the short range potential
    void Initialize(Graph* graph, Globaleventinfo* globevent); 
    double Calculate_longrange(int layer, bool cut_out_discs, myvec sim_box_size, Globaleventinfo* globevent); // Calculate long-range part of Coulomb interaction

    vector<double> layercharge;
    vector<double> longrange_cache;
    vector<double> positional_average;
  
    int number_of_layers;

private:
    vector< vector <double> > precalculate_disc_contrib; // Precalculated disc contributions
    double Calculate_disc_contrib(int calculate_layer, int contrib_layer, myvec sim_box_size, Globaleventinfo* globevent); // Calculate disc contributions
  
    vector<int> first_contributing_layer; // What is the first layer that contributes to the relevant layer?
    vector<int> final_contributing_layer; // What is the last layer that contributes to the relevant layer?
  
};

void Longrange::Update_cache(myvec sim_box_size, Globaleventinfo* globevent) {
    for (int i=0; i<number_of_layers; i++) {
        longrange_cache[i] = Calculate_longrange(i,true, sim_box_size, globevent);
    }
}

double Longrange::Get_cached_longrange(int layer) {
    return longrange_cache[layer];
}

void Longrange::Reset() {
    for (int i=0; i<number_of_layers; i++) {
        layercharge[i] = 0.0;
        longrange_cache[i] = 0.0;
    }
}

void Longrange::Initialize (Graph* graph, Globaleventinfo* globevent) {

    number_of_layers = ceil(graph->sim_box_size.x()/graph->hopdist);
 
    vector<double> positional_sum;
    vector<int> number_of_charges;
    vector<bool> flagged_for_deletion;
    
    for (int ilayer=0; ilayer<number_of_layers; ilayer++) {
        positional_sum.push_back(0.0);
        number_of_charges.push_back(0);
        flagged_for_deletion.push_back(false);
    }        
        
    for (unsigned int inode=0; inode<graph->nodes.size(); inode++) {
        double posx = graph->nodes[inode]->node_position.x();
        int iposx = floor(posx/graph->hopdist);
        positional_sum[iposx] += posx;
        number_of_charges[iposx]++;
        graph->nodes[inode]->layer_index = iposx;
    }
    
    int rem_layers = 0;
    
    for (int ilayer=0; ilayer<number_of_layers; ilayer++) {
        if(number_of_charges[ilayer] != 0) {
            positional_average.push_back(positional_sum[ilayer]/number_of_charges[ilayer]);
        }
        else {
            flagged_for_deletion[ilayer] = true;
            rem_layers++;
        }
    }
    
    number_of_layers -= rem_layers;
    
    layercharge.resize(number_of_layers);
    longrange_cache.resize(number_of_layers);

    precalculate_disc_contrib.resize(number_of_layers);
    first_contributing_layer.resize(number_of_layers);
    final_contributing_layer.resize(number_of_layers);             
  
    for (int ilayer=0;ilayer<number_of_layers;ilayer++) {
        // define for every layer, how many other layers are within the coulomb cut off radius from this layer
        double define_layer = positional_average[ilayer];
    
        int start_index = 0;
        bool startfound = false;
        while (!startfound) {
            double start_layer = positional_average[start_index];
            if((define_layer-start_layer)<=globevent->coulcut) {
                startfound = true;
                first_contributing_layer[ilayer]=start_index;
            }
            start_index++;
        }
    
        int final_index = number_of_layers-1;
        bool finalfound = false;
        while (!finalfound) {
            double final_layer = positional_average[final_index];
            if((final_layer-define_layer)<=globevent->coulcut) {
                finalfound = true;
                final_contributing_layer[ilayer]=final_index;
            }
            final_index--;
        }
    
        int number_of_contributing_layers = final_contributing_layer[ilayer]-first_contributing_layer[ilayer]+1;
        precalculate_disc_contrib[ilayer].resize(number_of_contributing_layers);
    }
  
    for(int i=0; i<number_of_layers; i++) {
        layercharge[i] = 0.0;
        longrange_cache[i] = 0.0;
        int first_layer = first_contributing_layer[i];
        int final_layer = final_contributing_layer[i];
        for (int j=first_layer; j<=final_layer; j++) {
            precalculate_disc_contrib[i][j-first_layer] = Calculate_disc_contrib(i,j,graph->sim_box_size,globevent);
        }
    }
}

double Longrange::Calculate_disc_contrib(int calculate_layer, int contrib_layer, myvec sim_box_size,Globaleventinfo* globevent) {
 
    double calcpos = positional_average[calculate_layer];
    double contribpos = positional_average[contrib_layer];
    double rdist = contribpos-calcpos;
  
    double contrib = globevent->coulcut-fabs(rdist); // Direct contribution (no image), factor 2 pi is missing, included in compute_longrange   
    double radiussqr = globevent->coulcut*globevent->coulcut-rdist*rdist; //radius of contributing disc
    
    for (long i=0; i<globevent->nr_of_lr_images; i++) {
   
        // Calculate contribution from images
        long dist1;
        long dist2;
        long sign;
        if (ldiv(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = -1;
            dist1 = i*sim_box_size.x() + 2*contribpos - rdist;
            dist2 = i*sim_box_size.x() + 2*sim_box_size.x() - 2*contribpos + rdist;
        }
        else { // odd generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = 1; 
            dist1 = (i+1)*sim_box_size.x() + rdist;
            dist2 = (i+1)*sim_box_size.x() - rdist;
        }
        double diag1 = sqrt(radiussqr+dist1*dist1);
        double diag2 = sqrt(radiussqr+dist2*dist2);
        contrib += sign*(diag1-fabs(dist1)+diag2-fabs(dist2));
    }
    return contrib;
}

double Longrange::Calculate_longrange(int layer, bool cut_out_discs, myvec sim_box_size, Globaleventinfo* globevent) {
    // Potential is expressed in multiples of e/(4*pi*epsilon) (with e the elementary charge>0)
    double plate_contrib1 = 0.0;
    double disc_contrib = 0.0;
    double PI = 3.14159265358979323846264338327950288419716939937510;
    for(int i=0; i<layer; i++) {
        double charge_i = 1.0*layercharge[i];
        double position_i = 1.0*positional_average[i];
        plate_contrib1 += position_i*charge_i; // potential of a charged plate between two electrodes
        double distance = layer-position_i;
        int first_layer = first_contributing_layer[layer];
        if (distance<=globevent->coulcut) {
            // Cut out short-range sphere
            disc_contrib -= charge_i*precalculate_disc_contrib[layer][i-first_layer];
        }
    }
    double plate_contrib2 = 0.0;
    for(int i=layer; i<number_of_layers; i++) {
        double charge_i = 1.0*layercharge[i];
        double rel_position_i = 1.0*(sim_box_size.x()-positional_average[i]);
        plate_contrib2 += rel_position_i*charge_i; // potential of a charged plate between two electrodes
        double distance = positional_average[i]-layer;
        int first_layer = first_contributing_layer[layer];
        if (distance<=globevent->coulcut) {
            // Cut out short-range sphere
            disc_contrib -= charge_i*precalculate_disc_contrib[layer][i-first_layer];
        }
    }
    if (!cut_out_discs) { disc_contrib = 0.0; }
    double layerpos = positional_average[layer];
    return 4*PI*(plate_contrib1*(1-layerpos/sim_box_size.x()) + plate_contrib2*(layerpos/sim_box_size.x()) + 0.5*disc_contrib)/(sim_box_size.y()*sim_box_size.z());
}

}}

#endif
