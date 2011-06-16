#ifndef ENERGYCORR_H
#define	ENERGYCORR_H

#include "qmpair.h"
#include "qmcalculator.h"
#include <votca/tools/histogramnew.h>
/**
    \brief Calculates reduced site energy correlation function

Callname : "energycorr"

To quantify the degree of correlation in site energies due to electrostatics, one can calculate the spatial correlation function of \f$E_i\f$ and \f$E_j\f$ at a distance \f$r_{ij}\f$
\f[C(r_{ij}) = \frac{  \langle \left( E_i-\langle E\rangle \right)
                   \left( E_j-\langle E\rangle \right)\rangle}
                   {\langle\left( E_i -\langle E\rangle \right)^2\rangle}\f]
where \f$\langle E\rangle\f$ is the average site energy. \f$C(r_{ij})\f$ is zero if \f$E_i\f$ and \f$E_j\f$ are uncorrelated and 1 if they are fully correlated. For a system of randomly oriented point dipoles, the correlation function decays as \f$1/r\f$ at large distances

Reference : J.Chem.Phys. 129, 034709 (2008)
*/

class EnergyCorr : public QMCalculator
{
public:
    EnergyCorr() {};
    ~EnergyCorr() {};

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

#endif	/* ENERGYCORR_H */

