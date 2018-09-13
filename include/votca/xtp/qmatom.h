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

#ifndef VOTCA_XTP_QMATOM_H
#define	VOTCA_XTP_QMATOM_H

#include <votca/tools/vec.h> 
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
    
   QMAtom (int index,std::string element, double x, double y, double z)
            :_index(index), _type( element), _nuccharge(0), _ecpcharge(0),_partialcharge(0.0)
            {_pos=tools::vec(x,y,z);}
            
   QMAtom (int index,std::string element,const tools::vec& pos)
            :_index(index), _type(element ),_nuccharge(0), _ecpcharge(0),_partialcharge(0.0)
            {_pos=pos;}
   
   
   QMAtom ()
            :_index(0), _type(""),_nuccharge(0), _ecpcharge(0),_partialcharge(0.0)
            {_pos=tools::vec(0.0);}
       
   const tools::vec & getPos() const {return _pos;}
   
  
   void setPos(tools::vec position){_pos=position;}

   const std::string & getType() const { return _type;}
   
  int getAtomID()const{ return _index;}
   
   int getNuccharge() const{ return _nuccharge-_ecpcharge;}
       
   void setPartialcharge(double q){_partialcharge=q;}
   const double & getPartialcharge() const { return _partialcharge;}

private:
    
   int _index;
   std::string _type;
   tools::vec _pos;// Bohr
   int _nuccharge;//nuc charge is set in aobasis fill and ecpfill
   int _ecpcharge;
   double _partialcharge;

   
 public: 
   
   void WriteToCpt(CptLoc parent){
       CheckpointWriter w(parent);

       w(_index, "index");
       w(_type, "type");
       w(_pos, "pos");
       w(_nuccharge, "nuccharge");
       w(_ecpcharge, "ecpcharge");
       w(_partialcharge, "partialcharge");
   }

   void ReadFromCpt(CptLoc parent){
       CheckpointReader r(parent);

       r(_index, "index");
       r(_type, "type");
       r(_pos, "pos");
       r(_nuccharge, "nuccharge");
       r(_ecpcharge, "ecpcharge");
       r(_partialcharge, "partialcharge");
   }
};
    
}}

#endif	// VOTCA_XTP_QMATOM_H 

