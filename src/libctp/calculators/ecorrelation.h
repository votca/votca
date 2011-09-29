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

#ifndef ECORRELATION_H
#define	ECORRELATION_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/histogramnew.h>

namespace votca { namespace ctp {

/**
    \brief Normalized site energy correlation function

The spatial correlations in site energies can be quantified using the correlation function (J. Chem. Phys. 129, 034709, 2008)
\f[C(r_{ij}) = \frac{  \langle \left( E_i-\langle E\rangle \right)
                   \left( E_j-\langle E\rangle \right)\rangle}
                   {\langle\left( E_i -\langle E\rangle \right)^2\rangle}\f]
where \f$E_i\f$ and \f$E_j\f$ are site energies, \f$ r_{ij} \f$ is the distance between the sites, \f$ \langle E\rangle \f$ is the average site energy. \f$ C(r_{ij}) \f$ is zero if \f$ E_i \f$ and \f$ E_j \f$ are uncorrelated and 1 if they are fully correlated. For a system of randomly oriented point dipoles, the correlation function decays as \f$ 1/r \f$ at large distances.

Callname: ecorrelation

*/
class Ecorrelation : public QMCalculator
{
public:
    Ecorrelation() {};
    ~Ecorrelation() {};

    const char *Description() { return "Normalized site energy correlation function"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

    /// A matching function to overload the standard one
    bool MyMatchingFunction(Bead *bead1, Bead *bead2, const vec & r, const double notused);
    
private:
    // histogram for recording energy correlation function
    HistogramNew _hist_corr;
    // histogram for recording counts (needed to normalize the first one)
    HistogramNew _hist_count;
    // histogram for recording correlation function squared (needed to calculate error bars)
    HistogramNew _hist_corr2;
    double _min, _max;
    double _mean_energy, _stdev_energy;
    int _nbins;
    /// name of the output file
    string _outfile;
};

}}

#endif	/* ECORRELATION_H */

