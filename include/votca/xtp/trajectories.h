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

#ifndef __VOTCA_KMC_TRAJECTORIES_H_
#define __VOTCA_KMC_TRAJECTORIES_H_

#include <votca/xtp/event.h>
#include <iostream>
#include <iomanip>


namespace votca { namespace xtp {
  
using namespace std;

class Trajectories {
    
public:
    
    void Init_trajectory(string filename);
    void Update_trajectory(Event* chosenevent);
    void Print_header();
    void Print_trajectory(double simtime);
 
 
private:
    
ofstream traj_stream;
char traj_file[100];
votca::tools::vec trajectory;
    
};

void Trajectories::Init_trajectory(string filename){
    strcpy(traj_file, filename.c_str());
    traj_stream.open(traj_file);
    this->Print_header();
    trajectory = votca::tools::vec(0.0,0.0,0.0);
}

void Trajectories::Print_header(){
    traj_stream << "time[s]" << "\t" << "x[nm]" << "\t" << "y[nm]" << "\t" << "z[nm]" << endl;
    traj_stream.flush();
}

void Trajectories::Update_trajectory(Event* chosenevent){
    votca::tools::vec traj_hop = chosenevent->link()->r12();
    trajectory = trajectory + traj_hop;
}

void Trajectories::Print_trajectory(double simtime){
    traj_stream << simtime << "\t";
    traj_stream << trajectory.x() << "\t";
    traj_stream << trajectory.y() << "\t";
    traj_stream << trajectory.z() << "\t";
    traj_stream << endl;
}
  
}}

#endif
