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

#ifndef __VOTCA_KMC_MESH_H_
#define __VOTCA_KMC_MESH_H_


namespace votca { namespace kmc {

template<class TObject>    
class Mesh {
public:

    /// Initialization of Mesh, which is nothing else than a vector (of variable dimensionality) of lists
    Mesh(double size_x, double size_y, votca::tools::vec simboxsize){ // Two-dimensional mesh of dimension (mesh_x x mesh_y)
        int dims_x = ceil(simboxsize.x()/size_x); // number of meshes
        int dims_y = ceil(simboxsize.y()/size_y);
        mesh_dims.clear(); mesh_dims.push_back(dims_x); mesh_dims.push_back(dims_y);
        mesh_sizes.clear(); mesh_sizes.push_back(size_x); mesh_sizes.push_back(size_y);
        mesh_vector.clear();
        std::cout << size_x << " " << size_y << " " << simboxsize << endl;
    };
    
    //Initialization of Mesh
    Mesh(double size_x, double size_y, double size_z, votca::tools::vec simboxsize){ // Three-dimensional mesh of dimension (mesh_x x mesh_y x mesh_z )
        int dims_x = ceil(simboxsize.x()/size_x); // number of meshes
        int dims_y = ceil(simboxsize.y()/size_y);
        int dims_z = ceil(simboxsize.z()/size_z);
        mesh_dims.clear(); mesh_dims.push_back(dims_x); mesh_dims.push_back(dims_y); mesh_dims.push_back(dims_z);
        mesh_sizes.clear(); mesh_sizes.push_back(size_x); mesh_sizes.push_back(size_y); mesh_sizes.push_back(size_z);
        mesh_vector.clear();
        std::cout << size_x << " " << size_y << " " << size_z << " " << simboxsize << endl;
    };
    
private:
    
    vector<TObject> mesh_vector;
    vector<int> mesh_dims;
    vector<double> mesh_sizes;
};

}} 

#endif

