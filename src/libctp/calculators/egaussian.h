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

#ifndef _EGAUSSIAN_H
#define	_EGAUSSIAN_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/average.h>

namespace votca { namespace ctp {

/** \brief Generates or modifies the Coulomb contribution to site energies.

Callname: egaussian

If a gaussian method is selected, site energies are generated from a Gaussian distribution of width sigma. If a correlation and a cutoff parameters are provided then a new site energy is computed for every site by averaging energies of the neighboring (up to a cutoff) sites. The new site energies are then uniformly scaled in order to reproduce the initial width (sigma) of the site energy distribution. If the shuffle method is selected, the (existing) site energies are randomly shuffled between the sites. This preserves the distribution of site energies but destroys spatial correlations.

*/

class Egaussian : public QMCalculator
{
public:
    Egaussian() {};
    ~Egaussian() {};

    const char *Description() { return "Generates or modifies Coulomb site energies."; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

    /// A matching function to overload the standard one
    bool MyMatchingFunction(Bead *bead1, Bead *bead2, const vec & r, const double notused);
    
private:
    /// sigma of the gaussian distribution
    double _sigma;
    /// true - correlated energies, false - uncorrelated energies
    bool _correl;
    /// correlation cutoff for generating correlated energies
    double _cutoff;
    /// temporary object to store energies
    std::map<CrgUnit *, Average<double> > _tmp_energy;

    bool (Egaussian::*_method)(QMTopology *top);

    bool AssignGaussian(QMTopology *top);
    bool AssignCorrelated(QMTopology *top);
    bool _shuffle(QMTopology *top);

};

}}

#endif	/* _EGAUSSIAN_H */

