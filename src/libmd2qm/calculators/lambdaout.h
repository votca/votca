/* 
 * File:   lambdaout.h
 * Author: mayfalk
 *
 * Created on May 20, 2011, 12:36 PM
 */

#ifndef _LAMBDAOUT_H
#define	_LAMBDAOUT_H




#include "qmpair.h"
#include "qmcalculator.h"


class CalcLambdaOut : public QMCalculator
{
public:
    CalcLambdaOut() {};
    ~CalcLambdaOut() {};

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
    Property * _options;
    // function pointer
    double (CalcLambdaOut::*_lambda_method)(const double &);

};


#endif	/* _LAMBDAOUT_H */

