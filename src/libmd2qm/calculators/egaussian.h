#ifndef _GENERATE_NRGS_H
#define	_GENERATE_NRGS_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/average.h>

/** \brief Gaussian-distributed site energies with or without spatial correlations

Callname: egaussian

Site energies of a gaussian distribution with width sigma are constructed. Then for every site a new energy is computed from the average of neighboring sites up to a cutoffradius. Then the new site energies are uniformly scaled in order to reproduce the intial sigma of the site energy distribution.

*/
class Egaussian : public QMCalculator
{
public:
    Egaussian() {};
    ~Egaussian() {};

    const char *Description() { return "Gaussian-distributed site energies with or without spatial correlations"; }

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


    void AssignGaussian(QMTopology *top);
    void AssignCorrelated(QMTopology *top);

};

#endif	/* _GENERATE_NRGS_H */

