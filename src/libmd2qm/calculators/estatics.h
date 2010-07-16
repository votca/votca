/* 
 * File:   Estatic_Calculator.h
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:23 AM
 */

#ifndef _CALC_ESTATICS_H
#define	_CALC_ESTATICS_H

#include "qmpair.h"
#include "qmcalculator.h"


class CalcEstatics : public QMCalculator
{
public:
    CalcEstatics() {};
    ~CalcEstatics() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);

    double CalcPot(Topology *atop, Molecule *mol);
    /// distance dependent dielectric constant
    double dist_dep_eps(const double &dist);
    /// constant epsilon
    double constant_epsilon(const double &dist);

private:
    // dielectric constant
    double _epsilon_dielectric;
    // parameter describing the decay of eps
    double _s_eps; 
    // function pointer 
    double (CalcEstatics::*_estatic_method)(const double &);
};

#endif	/* _CALC_ESTATICS_H */

