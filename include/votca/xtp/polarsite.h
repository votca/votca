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

#ifndef __VOTCA_XTP_POLARSITE_H
#define __VOTCA_XTP_POLARSITE_H


#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/segment.h>
namespace votca { namespace xtp {
    
class PolarSite
{

public:

    PolarSite(int id, const std::string& name, const Eigen::VectorXd& pos)
            : _id(id),              _name(name),_pos(pos),
            PhiP(0.0),PhiU(0.0){};
        
     
    // GET & SET & IMPORT FUNCTIONS
    int    getId() { return _id; }
    const std::string &getName() const{ return _name; }
    const Eigen::Vector3d   &getPos() const{ return _pos; }
    // Point charge 0,1 dipole,2 quad 
    int            getRank() { return _rank; }
   
  
    // COORDINATE TRANSFORMATION
    void            Translate(const Eigen::VectorXd& shift);
    void            Rotate(const  Eigen::Matrix3d &rot, const Eigen::VectorXd &refPos);
 
    Eigen::MatrixXd getQuad_cartesian();
    Eigen::VectorXd getQuad();//Q20, Q21c, Q21s, Q22c, Q22s
    Eigen::VectorXd getDipole();// x,y,z
    
    const Eigen::VectorXd & getMultipoles()const {return _multipoles;}//Q00,Q1x, Q1y, Q1z,Q20, Q21c, Q21s, Q22c, Q22s,
    
    void WriteMpsLine(std::ostream &out, string unit = "angstrom");

    void Induce(double wSOR);
    
    

private:

    int     _id;
    std::string  _name;
    Eigen::Vector3d _pos;
       
    int     _rank;

    Eigen::Matrix3d _Ps;
    Eigen::Vector3d _localpermanetField;
    Eigen::Vector3d _localinducedField;
    Eigen::Vector3d _inducedDipole;
    Eigen::Vector3d _inducedDipole_old;
    
    Eigen::VectorXd _multipoles; //Q00,Q1x, Q1y, Q1z,Q20, Q21c, Q21s, Q22c, Q22s,
    
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)

};


}}


#endif
