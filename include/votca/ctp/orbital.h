/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __VOTCA_CTP_ORBITAL_H
#define	__VOTCA_CTP_ORBITAL_H

namespace votca { namespace ctp {
using namespace votca::tools;

using namespace std;

/**
    \brief information about an orbital
 
    The Orbital class stores orbital id, energy, MO coefficients
    
*/
class Orbital 
{
public:   

    Orbital() { };
   ~Orbital() { };

    const int     &getId() const { return _id; }
    const double  &getEnergy() const { return _energy; }

    inline const double setEnergy( const double &energy ) { _energy = energy; }

protected:
    int                 _id;
    double              _energy;
    vector<double>      _coefficients;
        
};


}}

#endif	/* __VOTCA_CTP_ORBITAL_H */

