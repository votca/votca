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

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;

class Longrange {

public:
  void Update_cache(); // Update cached longrange contributions
  double Get_cached_longrange(int layer); // Return cached value
  void Reset();
  void Initialize(long number_of_long_range_images); //note that the number of images for the calculation of the long range potential should be considerably larger than the number for the short range potential
  double Calculate_longrange(int layer, bool cut_out_discs); // Calculate long-range part of Coulomb interaction

  vector<double> _layercharge;
  vector<double> _longrange_cache;
  vector<double> _positional_average;
  
  int _number_of_layers;
  int _number_of_longrange_images;
  double _device_length;
  double _coulomb_cut_off;
  double _box_dimension_Y;
  double _box_dimension_Z;

private:
  vector< vector <double>> precalculate_disc_contrib; // Precalculated disc contributions
  double Calculate_disc_contrib(int x, int rdist, long N); // Calculate disc contributions
  
  vector<int> first_contributing_layer; // What is the first layer that contributes to the relevant layer?
  vector<int> final_contributing_layer; // What is the last layer that contributes to the relevant layer?
  
  double pi;
};

void Longrange::Update_cache() {
  for (int i=0; i<_number_of_layers; i++) {
    _longrange_cache[i] = compute_longrange(i,true);
  }
}

double Longrange::Get_cached_longrange(int layer) {
  return _longrange_cache[layer];
}

void Longrange::Reset() {
  for (int i=0; i<_number_of_layers; i++) {
    _layercharge[i] = 0;
    _longrange_cache[i] = 0;
  }
}

void Longrange::Initialize (long number_of_long_range_images) {

  _layercharge.resize(_number_of_layers);
  _longrange_cache.resize(_number_of_layers);
  
  precalc_disc_contrib.resize(_number_of_layers);
  first_contributing_layer.resize(_number_of_layers);
  final_contributing_layer.resize(_number_of_layers);
  
  for (int ilayer=0;ilayer<_number_of_layers;ilayer++) {
    // define for every layer, how many other layers are within the coulomb cut off radius from this layer
    double define_layer = _positional_average[ilayer];
    
    int start_index = 0;
    bool startfound = false;
    while (!startfound) {
      double start_layer = _positional_average[start_index];
      if((define_layer-start_layer)<=_coulomb_cut_off) {
        startfound = true;
        first_contributing_layer[ilayer]=start_index;
      }
      start_index++;
    }
    
    int final_index = _number_of_layers-1;
    bool finalfound = false;
    while (!finalfound) {
      double final_layer = _positional_average[final_index];
      if((final_layer-define_layer)<=_coulomb_cut_off) {
        finalfound = true;
        final_contributing_layer[ilayer]=final_index;
      }
      final_index--;
    }
    
    int number_of_contributing_layers = final_contributing_layer[ilayer]-first_contributing_layer[ilayer]+1;
    precalc_disc_contrib[ilayer].resize(number_of_contributing_layers);
  }
   
  pi = 3.14159265358979323846264338327950288419716939937510;
  
  for(int i=0; i<_number_of_layers; i++) {
    _layercharge[i] = 0.0;
    _longrange_cache[i] = 0.0;
    int first_layer = first_contributing_layer[i];
    int final_layer = final_contributing_layer[i];
    for (int j=first_layer; j<=final_layer; j++) {
      precalculated_disc_contrib[i][j-first_layer)] = calculate_disc_contrib(i,j,number_of_long_range_images);
    }
  }
}

double Longrange::Calculate_disc_contrib(int calculate_layer, int contrib_layer, long number_of_long_range_images) {
  //direct contribution (no image), factor 2 pi is missing, included in compute_longrange   
  double calculate_position = _positional_average[calculate_layer];
  double contrib_position = _positional_average[contrib_layer];
  double rdist = contrib_position-calculate_position;
  
  double contrib = _coulomb_cut_off-abs(rdist); // Direct contribution
  double radiussqr = _coulomb_cut_off*_coulomb_cut_off-rdist*rdist; //radius of contributing disc
    
  for (long i=0; i<number_of_long_range_images; i++) {
  // Calculate contribution from images
    long dist1;
    long dist2;
    long sign;
    if (div(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
      sign = -1;
      dist1 = i*_device_length + 2*x - rdist;
      dist2 = i*_device_length + 2*_device_length - 2*x + rdist;
    }
    else { // odd generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
      sign = 1; 
      dist1 = (i+1)*_device_length + rdist;
      dist2 = (i+1)*_device_length - rdist;
    }
    double diag1 = sqrt(radiussqr+dist1*dist1);
    double diag2 = sqrt(radiussqr+dist2*dist2);
    contrib += sign*(diag1-abs(dist1)+diag2-abs(dist2));
  }
  return contrib;
}

double Longrange::Calculate_longrange(int layer, bool cut_out_discs) {
  // Potential is expressed in multiples of e/(4*pi*epsilon) (with e the elementary charge>0)
  double plate_contrib1 = 0.0;
  double disc_contrib = 0.0;
  for(int i=0; i<layer; i++) {
    double charge_i = 1.0*_layercharge[i];
    double position_i = 1.0*_positional_average[i];
    plate_contrib1 += position_i*charge_i; // potential of a charged plate between two electrodes
    double distance = layer-position_i;
    int first_layer = first_contributing_layer[layer];
    if (distance<=_coulomb_cut_off) {
      // Cut out short-range sphere
      disc_contrib -= charge_i*precalculated_disc_contrib[layer][i-first_layer];
    }
  }
  double plate_contrib2 = 0.0;
  for(int i=layer; i<_number_of_layers; i++) {
    double charge_i = 1.0*layercharge[i];
    double rel_position_i = 1.0*(_device_length-_positional_average[i]);
    plate_contrib2 += rel_position_i*charge_i; // potential of a charged plate between two electrodes
    double distance = _positional_average[i]-layer;
    int first_layer = first_contributing_layer[layer];
    if (distance<=_coulomb_cut_off) {
      // Cut out short-range sphere
      disc_contrib -= charge_i*precalc_disc_contrib[layer][i-first_layer];
    }
  }
  if (!cut_out_discs) { disc_contrib = 0.0; }
  double layerpos = _positional_average[layer];
  return 4*PI*(plate_contrib1*(1-layerpos/_device_length) + plate_contrib2*(layerpos/_device_length) + 0.5*disc_contrib)/(_box_dimension_Y*_box_dimension_Z);
}

#endif