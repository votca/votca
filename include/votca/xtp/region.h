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

#ifndef VOTCA_XTP_REGION_H
#define	VOTCA_XTP_REGION_H
#include <votca/xtp/eigen.h>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>

/**
* \brief base class to derive regions from
*
* 
* 
*/

namespace votca { namespace xtp {
   
template<class T>  class AtomContainer{
    public:
               
        AtomContainer(std::string name,int id):_name(name),_id(id),_position_valid(false){};
        
        const std::string& getName()const{return _name;}
        
        int getId()const{return _id;}
        
        int size()const{return _atomlist.size();}
        
        void push_back(const T& atom){_atomlist.push_back(atom);_position_valid=false;}

        const T& at(int index)const{return _atomlist.at(index);}
        T& at(int index){return _atomlist.at(index);}

        const T& operator[](int index)const{return _atomlist[index];}
        
        typename std::vector<T>::iterator begin(){return _atomlist.begin();}
        typename std::vector<T>::iterator end(){return _atomlist.end();}
        
        typename std::vector<T>::const_iterator begin()const{return _atomlist.begin();}
        typename std::vector<T>::const_iterator end()const{return _atomlist.end();}
        
        void calcPos(){
            tools::Elements element;
            _pos=Eigen::Vector3d;
            double totalmass=0.0;
            for (const T& atom:_atomlist){
                double mass=element.getMass(atom.getType())
                totalmass+=mass;
                _pos+=mass*atom.getPos()
            }
            _pos/=totalmass;
            _position_valid=true;
        }
       
        const Eigen::Vector3d& getPos()const{
            if(!_position_valid){calcPos();}
            return _pos;
        }
        
        void Translate(const Eigen::Vector3d& shift){
            for(const T& atom:_atomlist){
                atom.Translate(shift);
            }
            calcPos();
        }
        
        void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos){
            for(const T& atom:_atomlist){
                atom.Rotate(R,ref_pos);
            }
            calcPos();
        }
        
      virtual void WriteToCpt(CptLoc parent)const=0;
            
      virtual void ReadFromCpt(CptLoc parent)=0;
            
  protected:

    std::vector<T> _atomlist;
    std::string _name;
    int _id;
    
    private:
    mutable bool _position_valid;
    mutable Eigen::Vector3d _pos;
        
      };   
    
    
}}

#endif	// VOTCA_XTP_ATOMCONTAINER_H
