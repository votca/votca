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
    
protected:
    double _cutoff;
};

inline void Neighborlist::Initialize(QMTopology *top, Property *options){
    _cutoff = options->get("options.neighborlist.cutoff").as<double>();
}

inline bool Neighborlist::EvaluateFrame(QMTopology *top)
{
    top->nblist().Cleanup();
    top->nblist().setCutoff(_cutoff);
    BeadList list;
    list.Generate(*top, "*");
    top->nblist().Generate(list);
}

}}

#endif	/* __NEIGHBORLIST_H */

