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
#define __NEIGHBORLIST_H

#include <votca/ctp/qmcalculator.h>
#include <votca/tools/globals.h>
#include <votca/ctp/qmpair.h>

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

    // list of all sub-properties of options.neighborlist.segments
    list<Property *> segments = options->Select("options.neighborlist.segments");
    list<Property *>::iterator property_iterator;

    cout << "Neighbor list is based on the following cutoffs [type:type:cutoff]" << endl;

    // loop over these subproperites
    for (property_iterator = segments.begin(); property_iterator != segments.end(); ++property_iterator){
        string types = (*property_iterator)->get("type").as<string>();
        double  cutoff = (*property_iterator)->get("cutoff").as<double>();

        Tokenizer tok(types, " ");
        vector <string> segment_types;
        tok.ToVector(segment_types);

        // check if user provided two segment types
        if ( segment_types.size() != 2) std::runtime_error("error, two segment types separated by a space are needed for each cutoff");

        // sanity check is needed here to insure that segment type exists
        top->GetCrgUnitTypeByName(segment_types[0]);
        top->GetCrgUnitTypeByName(segment_types[1]);

        cout << " " << segment_types[0] << ":" << segment_types[1] << ":" << cutoff << "" << endl;

        _cutoff.push_back(cutoff);
        _typeA.push_back(segment_types[0]);
        _typeB.push_back(segment_types[1]);

    }
    
    cout <<"Using "<<_cutoff.size()<< " different types of pairs."<<endl;

  }

inline bool Neighborlist::EvaluateFrame(QMTopology *top)
{
    top->nblist().Cleanup();

    int pairtype;

    for (pairtype = 0; pairtype< _cutoff.size(); pairtype++) {

    //cout <<"TypeA:"<< _typeA[pairtype].c_str() << " TypeB:"<< _typeB[pairtype].c_str()<< " Cutoff:"<<_cutoff[pairtype]<<endl;
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
       if (tools::globals::verbose) {
       QMNBList& nblist=top->nblist();
        for (QMNBList::iterator ipair = nblist.begin(); ipair != nblist.end(); ++ipair) {

        QMPair *pair = *ipair;
        QMCrgUnit *crg1 = pair->Crg1PBCCopy();
        QMCrgUnit *crg2 = pair->Crg2PBCCopy();
        cout<<" id segment A:" << crg1->getId() <<"    id segment B:" << crg2->getId()<< "    distance:"<< pair->dist()<< endl;
        //cout<<"type A:" << crg1->getType().GetName() <<" name A:" << crg1->getName() << " id A:" << crg1->getId()<< endl;
    }
    }
}

}}

#endif  /* __NEIGHBORLIST_H */