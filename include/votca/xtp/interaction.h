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

#ifndef __VOTCA_XTP_INTERACTOR_H
#define __VOTCA_XTP_INTERACTOR_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/polarsite.h>
namespace votca { namespace xtp {\
    
  
    /**
    \brief Class to represent Atom/Site in electrostatic 

     The units are atomic units, e.g. Bohr, Hartree.
*/
class Interactor
{

public:
       
    double InteractStatic(StaticSite& Site1, StaticSite& Site2);
    
    double InteractInduction(PolarSite& Site1, PolarSite& Site2, double a=0.39);
        
 
protected:
       
    
    Eigen::MatrixXd FillTholeInteraction(const PolarSite& Site1, PolarSite& Site2, double a);
    Eigen::MatrixXd FillInteraction(const PolarSite& Site1, PolarSite& Site2);
    
    double _energy=0;
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)
    Eigen::Vector3d _inducedDipole;

};


}}


#endif
