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

#ifndef _EOUTERSPHERE2_H
#define	_EOUTERSPHERE2_H

#include <votca/ctp/topology.h>
#include <votca/ctp/qmcalculator2.h>

namespace votca { namespace ctp {

/**
    \brief Outer-sphere reorganization energy

Callname: eoutersphere

*/
class Eoutersphere2 : public QMCalculator2
{
public:
    Eoutersphere2() {};
   ~Eoutersphere2() {};

    void Initialize(Topology *top, Property *options);
    bool EvaluateFrame(Topology *top);

    void ConstLambda(Topology *top);
    void SpheresLambda(Topology *top);
    void DielectricLambda(Topology *top);

private:
    string      _method;
    double      _pekarFactor;
    double      _lambdaConstant;
    double      _lambdaCutoff;
    double (Eoutersphere2::*_lambda_method)(const double &);

};

}}

#endif	/* _EOUTERSPHERE_H */

