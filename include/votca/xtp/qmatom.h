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

#ifndef __VOTCA_XTP_QMATOM_H
#define	__VOTCA_XTP_QMATOM_H

#include <votca/tools/vec.h> 
#include <votca/xtp/aoshell.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>


namespace votca { namespace xtp {
    
/**
 *    \brief container for QM atoms 
 *
 *    Stores atom type, coordinates, charge
 */    
class QMAtom
{
     friend class AOBasis;
public:
    
   QMAtom (int _index,std::string _element, double _x, double _y, double _z)
            :index(_index), type( _element), nuccharge(0), ecpcharge(0),partialcharge(0.0)
            {pos=tools::vec(_x,_y,_z);}
            
   QMAtom (int _index,std::string _element, tools::vec _pos)
            :index(_index), type( _element ),nuccharge(0), ecpcharge(0),partialcharge(0.0)
            {pos=_pos;}
   
   
   QMAtom ()
            :index(0), type(""),nuccharge(0), ecpcharge(0),partialcharge(0.0)
            {pos=tools::vec(0.0);}
       
   const tools::vec & getPos() const {return pos;}
   
  
   void setPos(tools::vec position){pos=position;}

   const std::string & getType() const { return type;}
   
  int getAtomID(){ return index;}
   
   int getNuccharge() { return nuccharge-ecpcharge;}
       
   void setPartialcharge(double _q){partialcharge=_q;}
   const double & getPartialcharge() const { return partialcharge;}
   
  
   
private:
    
   int index;
   std::string type;
   tools::vec pos;// Bohr
   int nuccharge;//nuc charge is set in aobasis fill and ecpfill
   int ecpcharge;
   double partialcharge;

   
 public: 
   
   void WriteToCpt(CptLoc parent){
       CheckpointWriter w(parent);

       w(index, "index");
       w(type, "type");
       w(pos, "pos");
       w(nuccharge, "nuccharge");
       w(ecpcharge, "ecpcharge");
       w(partialcharge, "partialcharge");
   }

   void ReadFromCpt(CptLoc parent){
       CheckpointReader r(parent);

       r(index, "index");
       r(type, "type");
       r(pos, "pos");
       r(nuccharge, "nuccharge");
       r(ecpcharge, "ecpcharge");
       r(partialcharge, "partialcharge");
   }
};
    
}}

#endif	/* __VOTCA_XTP_QMATOM_H */

