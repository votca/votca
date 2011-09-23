#ifndef _ECOULOMB_H
#define	_ECOULOMB_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

/**
    \brief Coulomb contribution to site energies

Callname: ecoulomb

Electrostatic contribution \f$E^\textrm{el}_i\f$ to site energies is computed directly from the Coulomb's law. For all sites \f$i\f$
    \f[E^\textrm{el}_{i}  = \frac{1}{4 \pi \epsilon_0} \sum_{a_i} \sum_{b_k, k\neq i }
\frac{ \left( q^c_{a_i} - q^n_{a_i} \right) q^n_{b_k}}{ \epsilon_{s} r_{a_i b_k}}\f]
where \f$r_{a_i b_k}=|\vec{r}_{a_i} - \vec{r}_{b_k}|\f$
is the distance between atoms \f$a_i\f$ (which are in site \f$i\f$) and \f$b_k\f$ (which are in site \f$k\f$). The minimum image convention is used to account for periodic boundary conditions. \f$q^n\f$ and \f$q^c\f$ are partial charges of a neutral and charged hopping site, \f$\epsilon_s\f$ is the static relative dielectric constant. Depending on the method, \f$\epsilon_s\f$ is either constant or distance-dependent. In the latter case screening length \f$s\f$ must be specified.

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
    // function pointer 
    double (Ecoulomb::*_estatic_method)(const double &);
};
}}

#endif	/* _ECOULOMB_H */

