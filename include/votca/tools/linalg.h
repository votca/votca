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

    /**
     * \brief solves A*x=b
     * @param x storage for x
     * @param A symmetric positive definite matrix for linear system
     * @param b inhomogenity
     * @param if A is not sysmetrix positive definite throws error code GSL_EDOM
     *
     * This function wrapps the cholesky linear system solver
     */
    void linalg_cholesky_solve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b);
    /**
     * \brief solves A*x=b
     * @param x storage for x
     * @param A matrix for linear equation system
     * @param b inhomogenity
     * @param residual if non-zero, residual will be stored here
     *
     * This function wrapps the qrsolver
     */
    void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::vector<double> *residual=NULL);

    /**
     * \brief solves A*x=b under the constraint B*x = 0
     * @param x storage for x
     * @param A matrix for linear equation system
     * @param b inhomogenity
     * @param constr constrained condition B (or is it the transposed one? check that)
     *
     * This function wrapps the qrsolver under constraints
     */
    void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr);

}}


#endif	/* _LINSOLVE_H */

