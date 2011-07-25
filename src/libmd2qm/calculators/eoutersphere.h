/* 
 * File:   lambdaout.h
 * Author: mayfalk
 *
 * Created on May 20, 2011, 12:36 PM
 */

/**                                  
    \brief Outer-sphere reorganization energy         

Callname: eoutersphere 

Outer-sphere reorganization energy, \f$ \lambda_{ij}^\textrm{out} \f$ is be calculated from the electric displacement fields created by the charge transfer complex
\f[ \lambda_{ij}^\textrm{out}=
\frac{c_p}{2\epsilon_0}\int_{V^\textrm{out}}d V
\left[ \vec{D}_I(\vec{r}) - \vec{D}_F(\vec{r}) \right]^2 \f]
where \f$ \epsilon_0 \f$ is the the permittivity of free space, \f$ \vec{D}_{I,F}(\vec{r}) \f$ are the electric displacement fields created by the charge transfer complex in the initial (charge on molecule \f$i\f$) and final (charge transferred to molecule \f$ j \f$) states,  \f$ V^\textrm{out} \f$ is the volume outside the complex, and \f$ c_p=\frac{1}{\epsilon_\textrm{opt}}-\frac{1}{\epsilon_\textrm{s}} \f$ is the Pekar factor, which is determined by the low (\f$ \epsilon_\textrm{s} \f$) and high (\f$ \epsilon_\textrm{opt} \f$) frequency dielectric permittivities.

The difference of the displacement fields at the position of an atom \f$ b_k \f$ outside the charge transfer complex (molecule \f$ k \ne i,j \f$)  is expressed as
\f[ \vec{D}_I(\vec{r}_{b_k}) - \vec{D}_F(\vec{r}_{b_k})  
  =  \sum_{a_i} \frac{q_{a_i}^c - q_{a_i}^n}{4\pi } \frac{ (\vec{r}_{b_k} - \vec{r}_{a_i} ) }
  {|\vec{r}_{b_k}-\vec{r}_{a_i}|^3} + \sum_{a_j} \frac{q_{a_j}^n - q_{a_j}^c}{4\pi } \frac{ (\vec{r}_{b_k}-\vec{r}_{a_j} ) }
  {|\vec{r}_{b_k}-\vec{r}_{a_j}|^3}
\f]
where \f$ q^n_{a_i} \f$ (\f$ q^c_{a_i} \f$) is the partial charge of atom \f$ a \f$ of the neutral (charged) molecule \f$ i \f$ in vacuum. The partial charges of neutral and charged molecules are obtained by fitting their values to reproduce the electrostatic potential of a single molecule (charged or neutral) in vacuum.

*/                                         

#ifndef _LAMBDAOUT_H
#define	_LAMBDAOUT_H




#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>


class CalcLambdaOut : public QMCalculator
{
public:
    CalcLambdaOut() {};
    ~CalcLambdaOut() {};

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
    double (CalcLambdaOut::*_lambda_method)(const double &);

};


#endif	/* _LAMBDAOUT_H */

