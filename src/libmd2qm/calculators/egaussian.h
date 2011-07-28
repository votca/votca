#ifndef _EGAUSSIAN_H
#define	_EGAUSSIAN_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/average.h>

/** \brief Generates or modifies the Coulomb contribution to site energies.

Callname: egaussian

If a gaussian method is selected, site energies are generated from a Gaussian distribution of width sigma. If a correlation and a cutoff parameters are provided then a new site energy is computed for every site by averaging energies of the neighboring, up to a cutoff, sites. The new site energies are then uniformly scaled in order to reproduce the initial width (sigma) of the site energy distribution. If the shuffle method is selected, the (existing) site energies are randomly shuffled between the sites. This preserves the distribution of site energies but destroys spatial correlations.

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

#endif	/* _EGAUSSIAN_H */

