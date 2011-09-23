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

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmtopology.h>

namespace votca { namespace ctp {

QMPair::QMPair(QMCrgUnit *crg1, QMCrgUnit *crg2, QMTopology * top)
  : std::pair<QMCrgUnit *, QMCrgUnit *>(crg1,crg2), _rate_12(0.),_rate_21(0.), _in_database(false)
{
    vec crg1nm,crg2nm;
    crg1nm =  crg1->GetCom();
    crg2nm =  crg2->GetCom();
    _r = top->BCShortestConnection(crg1nm, crg2nm);
    _crg2 = second;

    // check if PBC:
    vec d = crg2nm - crg1nm;
    if (abs(d - _r) > 1e-8) {
        _ghost = new QMCrgUnit();
	_ghost->copyCrgUnit(*crg2);
        vec displ = (_r - d);
        _ghost->shift(displ);
        _crg2 = _ghost;
    }
    else {
        _ghost=NULL;
    }
}

double QMPair::calcJeff2(){
    vector<double>::iterator itj=_Js.begin();
    double j=0.;
    for (;itj!= _Js.end(); itj++){
        j+=(*itj)*(*itj);
    }
    j/= double(_Js.size());
    return j;
}

}}