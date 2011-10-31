/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __NEIGHBORLIST_H
#define	__NEIGHBORLIST_H

#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

/** \brief Constructs a list of neighboring conjugated segments.

Callname: neighborlist

Two segments are added to this list if the distance between centers of mass of any their rigid fragments is below a certain cutoﬀ. This allows neighbors to be selected on a criterion of minimum distance of approach rather than center of mass distance, which is useful for molecules with anisotropic shapes.

*/
class Neighborlist : public QMCalculator
{
public:
    Neighborlist() {}
    ~Neighborlist() {}

    const char *Description() { return "Constructs a list of neighboring conjugated segments"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    
//protected:
private:
    vector<double> _cutoff;
    vector<string> _typeA;
    vector<string> _typeB;

};

inline void Neighborlist::Initialize(QMTopology *top, Property *options){
    //DENIS ??? _cutoff = options->get("options.neighborlist.cutoff").as<double>();
    _cutoff.push_back(0.7);
    _typeA.push_back("B1295");
    _typeB.push_back("B1295");
    
}

inline bool Neighborlist::EvaluateFrame(QMTopology *top)
{
    top->nblist().Cleanup();

    int pairtype;

    cout <<"There are "<<_cutoff.size()<< " different types of pairs for the neighborlist cutoff"<<endl;

    for (pairtype = 0; pairtype< _cutoff.size(); pairtype++) {
    cout <<"TypeA:"<< _typeA[pairtype].c_str() << " TypeB:"<< _typeB[pairtype].c_str()<< " Cutoff:"<<_cutoff[pairtype]<<endl;
    top->nblist().setCutoff(_cutoff[pairtype]);
    BeadList list1, list2;
    list1.Generate(*top, _typeA[pairtype]);
    list2.Generate(*top, _typeB[pairtype]);
     
      // is it same types or different types?
        if(_typeA[pairtype]== _typeB[pairtype])
            top->nblist().Generate(list1);
        else
        top->nblist().Generate(list1, list2);
    }
}

}}

#endif	/* __NEIGHBORLIST_H */

