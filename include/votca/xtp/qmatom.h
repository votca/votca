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

#ifndef __VOTCA_XTP_QMATOM_H
#define	__VOTCA_XTP_QMATOM_H

// Binary archive 
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
    
namespace votca { namespace xtp {
    
/**
 *    \brief container for QM atoms 
 *
 *    Stores atom type, coordinates, charge
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
    
}}

#endif	/* __VOTCA_XTP_QMATOM_H */

