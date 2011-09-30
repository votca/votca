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

#ifndef _ECOULOMB_H
#define	_ECOULOMB_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

/**
    \brief Electrostatic contribution to site energies

Callname: ecoulomb

*/
class Ecoulomb : public QMCalculator
{
public:
    Ecoulomb() {};
    ~Ecoulomb() {};

    const char *Description() { return "Coulomb contribution to site energies"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

    double CalcPot(Topology *atop, Molecule *mol);
    /// distance dependent dielectric constant
    double dist_dep_eps(const double &dist);
    /// constant dielectric constant
    double constant_epsilon(const double &dist);

private:
    // dielectric constant
    double _epsilon_dielectric;
    // screening length of the dielctirc constant
    double _s_eps;
    Property * _options;
    //if you want to compute estatics using a cutoff
    bool _has_cutoff;
    //cutoff in nm
    double _ecoulomb_cutoff;

    // function pointer 
    double (Ecoulomb::*_estatic_method)(const double &);
};
}}

#endif	/* _ECOULOMB_H */

