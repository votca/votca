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

#ifndef VOTCA_XTP_ATOMCONTAINER_H
#define	VOTCA_XTP_ATOMCONTAINER_H
#include <votca/xtp/eigen.h>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>
#include <limits>

/**
* \brief Basic Container for QMAtoms,PolarSites and Atoms
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
        T& operator[](int index){return _atomlist[index];}
        
        typename std::vector<T>::iterator begin(){return _atomlist.begin();}
        typename std::vector<T>::iterator end(){return _atomlist.end();}
        
        typename std::vector<T>::const_iterator begin()const{return _atomlist.begin();}
        typename std::vector<T>::const_iterator end()const{return _atomlist.end();}
        
        
       
        const Eigen::Vector3d& getPos()const{
            if(!_position_valid){calcPos();}
            return _pos;
        }

        //calculates the lowest and highest point in the cube, sorrounding the molecule
        std::pair<Eigen::Vector3d,Eigen::Vector3d> CalcSpatialMinMax() const{
            std::pair<Eigen::Vector3d,Eigen::Vector3d> result;
            Eigen::Vector3d min=std::numeric_limits<double>::max()*Eigen::Vector3d::Ones();
            Eigen::Vector3d max=std::numeric_limits<double>::min()*Eigen::Vector3d::Ones();
            for (const T& atom : _atomlist){
                const Eigen::Vector3d& pos=atom.getPos();
                if (pos.x()<min.x()) min.x()=pos.x();
                if (pos.x()>max.x()) max.x()=pos.x();
                if (pos.y()<min.y()) min.y()=pos.y();
                if (pos.y()>max.y()) max.y()=pos.y();
                if (pos.z()<min.z()) min.z()=pos.z();
                if (pos.z()>max.z()) max.z()=pos.z();
            }
            result.first=min;
            result.second=max;
            return result;
        }

        std::vector<std::string> FindUniqueElements()const{
            std::vector<std::string> result;
            for (const T& atom : _atomlist){
                if(std::find(result.begin(), result.end(), atom.getElement()) == result.end()) {
                    result.push_back(atom.getElement());
                }
            }
            return result;
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
        
    void WriteToCpt(CheckpointWriter& w)const{
        w(_name,"name");
        w(_id,"id");
    for (unsigned i=0;i<_atomlist.size();i++) {
        CheckpointWriter s = w.openChild( T::Identify() + std::to_string(i));
        _atomlist[i].WriteToCpt(s);
    }    
    }

    void ReadFromCpt(CheckpointReader& r){
        r(_name,"name");
        r(_id,"id");
    size_t count = r.getNumDataSets();
    _atomlist.resize(0);
    _atomlist.reserve(count);
    for (size_t i = 0; i < count; ++i) {
       CheckpointReader c = r.openChild( T::Identify() + std::to_string(i));
        T temp;
        temp.ReadFromCpt(c);
        _atomlist.emplace_back(temp);
    }      
    }
            
protected:

    std::vector<T> _atomlist;
    std::string _name;
    int _id;
    
private:
    mutable bool _position_valid;
    mutable Eigen::Vector3d _pos;

    void calcPos() const{
        tools::Elements element;
        _pos=Eigen::Vector3d::Zero();
        double totalmass=0.0;
        for (const T& atom:_atomlist){
            double mass=element.getMass(atom.getElement());
            totalmass+=mass;
            _pos+=mass*atom.getPos();
        }
        _pos/=totalmass;
        _position_valid=true;
    }
        
      };   
    
    
}}

#endif	// VOTCA_XTP_ATOMCONTAINER_H
