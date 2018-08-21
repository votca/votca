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

#ifndef _VOTCA_XTP_QMSTATE_H
#define _VOTCA_XTP_QMSTATE_H

#include<string>


namespace votca {
namespace xtp {
    
    
    class QMStateType{
    public:
        enum statetype {Singlet,Triplet,PQPstate,DQPstate,KSstate,Gstate};
        
        statetype Type()const{return _type;}
        
        void FromString(const std::string& statetypestring);
        
        std::string ToString();
    
        std::string ToLongString();
        
        bool operator ==(const QMStateType& rhs) const{
            return _type==rhs.Type();
        }
        
        bool operator==(const QMStateType::statetype& rhs){
            return _type==rhs;
        }
    
    private:
        
        statetype _type;               
    };
/**
 *  \brief  Identifier for QMstates. Strings like S1 are converted into enum +zero indexed int
 *
 *
 */
   
class QMState {

public:
    void FromString(const std::string& statestring);
    
    std::string ToString();
    
    std::string ToLongString();
    
    const QMStateType& Type()const{return _type;}
    
    bool isTransition()const{return _transition;}
    int Index()const{return _index;}
    
    bool operator ==(const QMState& rhs) const{
        return (_type==rhs.Type() && _index==rhs.Index());
        }
    
private:
 
    QMStateType _type;
    
    int _index;
    
    bool _transition;
    
    

};

}
}

#endif /* _VOTCA_XTP_QMSTATE_H */
