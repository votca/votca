/* 
 * File:   linsolve.h
 * Author: ruehle
 *
 * Created on March 15, 2010, 12:09 PM
 */

#ifndef __VOTCA_LINSOLVE_H
#define	__VOTCA_LINSOLVE_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


namespace votca { namespace tools {
    namespace ub = boost::numeric::ublas;

    void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::vector<double> *residual=NULL);

    void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr);

}}


#endif	/* _LINSOLVE_H */

