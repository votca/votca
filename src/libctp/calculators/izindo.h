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

#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>

namespace votca { namespace ctp {

/**
	\brief Semi-empirical electronic coupling elements for all neighbor list pairs

Semi-emprirical (ZINDO) electronic coupling elements for all conjuageted segments from the neighbout list. Requires molecular orbitals in GAUSSIAN format.

Callname: izindo

*/
class Izindo : public PairCalculator
{
public:
    Izindo() {};
    ~Izindo() {};

    //const char *Description() { return "Semi-empirical electronic coupling elements for all neighbor list pairs"; }

    void EvaluatePair(QMTopology *top, QMPair *pair);
};

inline void Izindo::EvaluatePair(QMTopology *top, QMPair *pair){
    CrgUnit *crg1 = pair->Crg1PBCCopy();
    CrgUnit *crg2 = pair->Crg2PBCCopy();
    vector <double> Js = top->GetJCalc().CalcJ(*crg1, *crg2);
    pair->setJs(Js);
}

}}

#endif	/* _CALC_INTEGRALS_H */

