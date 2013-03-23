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

#ifndef __VOTCA_CTP_ORBITALS_H
#define	__VOTCA_CTP_ORBITALS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/version.hpp>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
/**
    \brief container for molecular orbitals
 
    The Orbital class stores orbital id, energy, MO coefficients
    
*/
class Orbitals 
{
public:   

    Orbitals();
   ~Orbitals();

    const int     &getBasisSetSize() const;
    const int     &getNumberOfLevels() const;
    const int     &getNumberOfElectrons() const;
        
    ub::symmetric_matrix<double>* getOverlap() { return &_overlap; }
    ub::matrix<double>* getOrbitals() { return &_mo_coefficients; }
    ub::vector<double>* getEnergies() { return &_mo_energies; }
    
    std::vector<int>* getDegeneracy( int level, double _energy_difference );
    
    bool ReadOrbitalsGaussian( const char * filename );
    bool ReadOverlapGaussian( const char * filename );
    bool ParseGaussianLog( const char * filename );
    bool Save( const char * filename );

    
protected:
    
    static const double                 _conv_Hrt_eV = 27.21138386;
    
    int                                 _basis_set_size;
    int                                 _occupied_levels;
    int                                 _unoccupied_levels;
    int                                 _electrons;
    
    bool                                _has_basis_set_size;
    bool                                _has_occupied_levels;
    bool                                _has_unoccupied_levels;
    bool                                _has_electrons;
    bool                                _has_degeneracy;
    
    std::vector<int>                    _active_levels;
    std::map<int, std::vector<int> >    _level_degeneracy;
     
    ub::vector<double>                  _mo_energies;    
    ub::matrix<double>                  _mo_coefficients;
    ub::symmetric_matrix<double>        _overlap;

private:

    /**
    * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
    * @return A map with key as a level and a vector which is a list of close lying orbitals
    */    
    bool CheckDegeneracy( double _energy_difference );
    
    // Allow serialization to access non-public data members
    friend class boost::serialization::access;
    
    // serialization itself
    template<typename Archive> 
    void serialize(Archive& ar, const unsigned version);
    
    // 

};

//BOOST_CLASS_VERSION(Orbitals, 1)
        
}}

#endif	/* __VOTCA_CTP_ORBITALS_H */

