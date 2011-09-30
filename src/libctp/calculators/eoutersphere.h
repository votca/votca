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

#ifndef _EOUTERSPHERE_H
#define	_EOUTERSPHERE_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

/**
    \brief Outer-sphere reorganization energy

Callname: eoutersphere

*/
class Eoutersphere : public QMCalculator
{
public:
    Eoutersphere() {};
    ~Eoutersphere() {};

    const char *Description() { return "Outer-sphere reorganization energy"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

    double CalcLambdaDielectric(Topology *atop, Molecule *mol, Bead *bead);
    double CalcLambdaDist(Topology *atop, QMPair *pair);
        /// constant epsilon
    void const_lambda(QMTopology *top);
    void spheres_lambda(QMTopology *top);
    void dielectric_lambda(QMTopology *top);

private:
    double _pekar;
    double _R_one, _R_two;
    double _lambda_const;
    double _lambda_cutoff;
    Property * _options;
    // function pointer
    double (Eoutersphere::*_lambda_method)(const double &);

};

}}

#endif	/* _EOUTERSPHERE_H */

