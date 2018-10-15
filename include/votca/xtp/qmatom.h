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

#include <votca/tools/elements.h>
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
    
   
            
   QMAtom (int index,std::string element,const Eigen::Vector3d pos)
            :_index(index), _element(element ),_pos(pos),_nuccharge(0)
            , _ecpcharge(0)
            {
                tools::Elements elements;
                _nuccharge=elements.getNucCrg(element);
            }

  QMAtom (int index,std::string element, double x, double y, double z):
            QMAtom(index,element,Eigen::Vector3d{x,y,z}){};
       
   
       
   const Eigen::Vector3d& getPos() const {return _pos;}
   
   void Translate(const Eigen::Vector3d& shift){
       _pos+=shift;
   }
   
    void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d &refPos) {
      Translate(-refPos);
      _pos = R*_pos;
      Translate(refPos); //Rotated Position
    }
   
   void setPos(const Eigen::Vector3d& position){_pos=position;}

   const std::string & getElement() const { return _element;}
   
  int getAtomID()const{ return _index;}
   
   int getNuccharge() const{ return _nuccharge-_ecpcharge;}

private:
    
   int _index;
   std::string _element;
   Eigen::Vector3d _pos;// Bohr
   int _nuccharge;//nuc charge is set in aobasis fill and ecpfill
   int _ecpcharge;


 public:

   void WriteToCpt(CheckpointWriter& w){
       w(_index, "index");
       w(_element, "element");
       w(_pos, "pos");
       w(_nuccharge, "nuccharge");
       w(_ecpcharge, "ecpcharge");
   }

   void ReadFromCpt(CheckpointReader& r){
       r(_index, "index");
       r(_element, "element");
       r(_pos, "pos");
       r(_nuccharge, "nuccharge");
       r(_ecpcharge, "ecpcharge");
   }
};
    
}}

#endif	// VOTCA_XTP_QMATOM_H 

