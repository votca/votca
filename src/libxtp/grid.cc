/* 
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/grid.h>
#include <math.h>       /* ceil */
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>


namespace votca { namespace xtp {

void Grid::printGridtoxyzfile(std::string filename){
        //unit is Angstrom in xyz file 
        std::ofstream points;
        points.open(filename.c_str(), std::ofstream::out);
        points << _gridpoints.size() << std::endl;
        points << std::endl;
        for ( const auto& point:_gridpoints){
            points << "X " << point.x()*tools::conv::bohr2ang << " "
                    << point.y()*tools::conv::bohr2ang << " "
                    << point.z()*tools::conv::bohr2ang << std::endl;
        }
        points.close();
        return;
    }    


void Grid::setupgrid(const QMMolecule& Atomlist){
            
    tools::Elements elements;
    std::pair<Eigen::Vector3d,Eigen::Vector3d> extension=Atomlist.CalcSpatialMinMax();
    Eigen::Array3d min=extension.first.array();
    Eigen::Array3d max=extension.second.array();
    Eigen::Array3d doublesteps=(max-min+2*_padding)/_gridspacing;
    Eigen::Array3i steps=(doublesteps.ceil()).cast<int>();

    // needed to symmetrize grid around molecule
    Eigen::Array3d padding=(doublesteps-steps.cast<double>())*_gridspacing*0.5+_padding;
    Eigen::Array3d minpos=min-padding;
    for(int i=0;i<=steps.x();i++){
        double x=minpos.x()+i*_gridspacing;
        for(int j=0;j<=steps.y();j++){
            double y=minpos.y()+j*_gridspacing;
            for(int k=0;k<=steps.z();k++){
                double z=minpos.z()+k*_gridspacing;
                bool is_valid = false;
                Eigen::Vector3d gridpos(x,y,z);
                    for (const QMAtom& atom : Atomlist){
                        const Eigen::Vector3d& atompos=atom.getPos();
                        double distance2=(gridpos-atompos).squaredNorm();
                        double atomcutoff=elements.getVdWChelpG(atom.getElement())*tools::conv::ang2bohr;
                        if ( distance2<(atomcutoff*atomcutoff)){
                            is_valid = false;
                            break;
                            }
                        else if ( distance2<(_cutoff*_cutoff))  is_valid = true;
                    }
                    if (is_valid){
                      _gridpoints.push_back(gridpos);
                }                          
            }                  
        }
    }

    _gridvalues=Eigen::VectorXd::Zero(_gridpoints.size());
    return;
}
    
    
}}
