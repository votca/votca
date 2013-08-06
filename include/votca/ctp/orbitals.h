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
#include <boost/numeric/ublas/io.hpp>
#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

// Text archive that defines boost::archive::text_oarchive
// and boost::archive::text_iarchive
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// XML archive that defines boost::archive::xml_oarchive
// and boost::archive::xml_iarchive
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

// XML archive which uses wide characters (use for UTF-8 output ),
// defines boost::archive::xml_woarchive
// and boost::archive::xml_wiarchive
#include <boost/archive/xml_woarchive.hpp>
#include <boost/archive/xml_wiarchive.hpp>

// Binary archive that defines boost::archive::binary_oarchive
// and boost::archive::binary_iarchive
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
/**
    \brief container for basic atoms 
     Stores atom type, coordinates, charge
 */    
class QMAtom
{
public:
    
   QMAtom (std::string _type, double _x, double _y, double _z, double _charge, bool _from_environment)
            : type( _type ), x(_x), y(_y), z(_z), charge(_charge), from_environment( _from_environment )
            {};
            
    QMAtom ()
            : type( "" ), x(0), y(0), z(0), charge(0), from_environment( false )
            {};     
            
            
   
            
            
   std::string type;
   double x;
   double y;
   double z;
   double charge;
   bool   from_environment;
   
   template<typename Archive> 
   void serialize(Archive& ar, const unsigned version) {
       ar & type;
       ar & x;
       ar & y;
       ar & z;
       ar & charge;
       ar & from_environment;
   }  
};
    
/**
    \brief container for molecular orbitals
 
    The Orbitals class stores orbital id, energy, MO coefficients
    
*/
class Orbitals 
{
public:   

    Orbitals();
   ~Orbitals();

    bool           hasBasisSetSize() { return _has_basis_set_size; }
    int            getBasisSetSize() { return (_has_basis_set_size) ? _basis_set_size : 0; }
    void           setBasisSetSize( const int &basis_set_size );
    
    int            getNumberOfLevels() { return (_has_occupied_levels && _has_unoccupied_levels) ? ( _occupied_levels + _unoccupied_levels ) : 0; }
    void           setNumberOfLevels( const int &occupied_levels, const int &unoccupied_levels );
    
    bool           hasNumberOfElectrons() { return _has_number_of_electrons; }
    int            getNumberOfElectrons() { return (_has_number_of_electrons) ? _number_of_electrons : 0; } ;
    void           setNumberOfElectrons( const int &electrons );
    
    ub::symmetric_matrix<double>* getOverlap() { return &_overlap; }
    ub::matrix<double>* getOrbitals() { return &_mo_coefficients; }
    ub::vector<double>* getEnergies() { return &_mo_energies; }

    ub::matrix<double>* getIntegrals() { return _integrals; }
    void setIntegrals( ub::matrix<double>* integrals ) { _has_integrals = true;  _integrals = integrals; }
 
    double getEnergy( int level) { return (_has_mo_energies) ? _conv_Hrt_eV*_mo_energies[level-1] : 0; }
    
    std::vector<int>* getDegeneracy( int level, double _energy_difference );
    std::vector< QMAtom* >* getAtoms() { return &_atoms; }
    
    bool hasSelfEnergy() { return _has_self_energy; }
    double getSelfEnergy() { return (_has_self_energy) ? _self_energy : 0; }

    bool hasQMEnergy() { return _has_qm_energy; }
    double getQMEnergy() { return (_has_qm_energy) ? _qm_energy : 0; }
    
    // returns indeces of a re-sorted in a descending order vector of energies
    void SortEnergies( std::vector<int>* index );
    
    QMAtom* AddAtom (std::string _type, double _x, double _y, double _z, double _charge = 0, bool _from_environment = false){
        //std::cout << _type << std::endl;
        QMAtom* pAtom = new QMAtom(_type, _x, _y, _z, _charge, _from_environment);
        _atoms.push_back( pAtom );
        return pAtom;
    }
    
    void setStorage( bool _store_orbitals, bool _store_overlap,  bool _store_integrals ) {
        _has_mo_coefficients = _store_orbitals;
        _has_overlap = _store_overlap;
        _has_integrals = _store_integrals;
    } 
        

private:
    
    static const double                 _conv_Hrt_eV = 27.21138386;

    bool                                _has_basis_set_size;
    int                                     _basis_set_size;   
    
    bool                                _has_occupied_levels;
    int                                     _occupied_levels;
    
    bool                                _has_unoccupied_levels;
    int                                     _unoccupied_levels;
    
    bool                                _has_number_of_electrons;
    int                                     _number_of_electrons;
    
    bool                                _has_level_degeneracy;
    std::map<int, std::vector<int> >        _level_degeneracy;
    
    bool                                _has_mo_energies;
    ub::vector<double>                      _mo_energies; 
    
    bool                                _has_mo_coefficients;
    ub::matrix<double>                      _mo_coefficients;
    
    bool                                _has_overlap;
    ub::symmetric_matrix<double>            _overlap;
    
    bool                                _has_charges;
    bool                                _has_atoms;
    std::vector< QMAtom* >                  _atoms;   

    bool                                _has_qm_energy;
    double                                  _qm_energy;
    
    bool                                _has_self_energy;
    double                                  _self_energy;
    
    bool                                _has_integrals;
    ub::matrix<double>*                 _integrals;

private:

    /**
    * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
    * @return A map with key as a level and a vector which is a list of close lying orbitals
    */    
    bool CheckDegeneracy( double _energy_difference );
    
    // Allow serialization to access non-public data members
    friend class boost::serialization::access;
    
    //Allow Gaussian object to access non-public data members
    friend class Gaussian;
    friend class Turbomole;
    
    // serialization itself (template implementation stays in the header)
    template<typename Archive> 
    void serialize(Archive& ar, const unsigned version) {

       ar & _has_basis_set_size;
       ar & _has_occupied_levels;
       ar & _has_unoccupied_levels;
       ar & _has_number_of_electrons;
       ar & _has_level_degeneracy;
       ar & _has_mo_energies;
       ar & _has_mo_coefficients;
       ar & _has_overlap;
       ar & _has_atoms;
       ar & _has_qm_energy;
       ar & _has_self_energy;
       ar & _has_integrals;
       
       if ( _has_basis_set_size ) { ar & _basis_set_size; }
       if ( _has_occupied_levels ) { ar & _occupied_levels; }
       if ( _has_unoccupied_levels ) { ar & _unoccupied_levels; }
       if ( _has_number_of_electrons ) { ar & _number_of_electrons; }
       if ( _has_level_degeneracy ) { ar & _level_degeneracy; }
       if ( _has_mo_energies ) { ar & _mo_energies; }
       if ( _has_mo_coefficients ) { ar & _mo_coefficients; }
       if ( _has_overlap ) { 
           // symmetric matrix does not serialize by default
            if (Archive::is_saving::value) {
                unsigned size = _overlap.size1();
                ar & size;
             }

            // copy the values back if loading
            if (Archive::is_loading::value) {
                unsigned size;
                ar & size;
                _overlap.resize(size);
             }
            
           for (unsigned i = 0; i < _overlap.size1(); ++i)
                for (unsigned j = 0; j <= i; ++j)
                    ar & _overlap(i, j); 
       }
       
       if ( _has_atoms ) { ar & _atoms; }
       if ( _has_qm_energy ) { ar & _qm_energy; }
       if ( _has_self_energy ) { ar & _self_energy; }     
       if ( _has_integrals ) { ar & _integrals; } 
    }
    
};

//BOOST_CLASS_VERSION(Orbitals, 1)
        
}}

#endif	/* __VOTCA_CTP_ORBITALS_H */

