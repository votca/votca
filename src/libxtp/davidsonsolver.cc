/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <stdexcept>
#include <iostream>

#include <votca/xtp/eigen.h>
#include <votca/xtp/davidsonsolver.h>



using boost::format;
using std::flush;

namespace votca { 
    namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver(ctp::Logger &log) : _log(log) { }

void DavidsonSolver::set_correction(std::string method) {
    if (method == "DPR") this->davidson_correction = CORR::DPR;
    else if (method == "OLSEN") this->davidson_correction = CORR::OLSEN;
    else throw std::runtime_error(method + " is not a valid Davidson correction method");
}

Eigen::ArrayXd DavidsonSolver::_sort_index(Eigen::VectorXd &V) const
{
    Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(V.rows(),0,V.rows()-1);
    std::sort(idx.data(),idx.data()+idx.size(),
              [&](int i1, int i2){return V[i1]<V[i2];});
    return idx; 
}

Eigen::MatrixXd DavidsonSolver::_get_initial_eigenvectors(Eigen::VectorXd &d, int size_initial_guess) const
{
    Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(d.size(),size_initial_guess);
    Eigen::ArrayXd idx = DavidsonSolver::_sort_index(d);

    for (int j=0; j<size_initial_guess;j++) {
        guess(idx(j),j) = 1.0; 
    }

    return guess;
}

Eigen::VectorXd DavidsonSolver::_dpr_correction(Eigen::VectorXd &r, Eigen::VectorXd &D, double lambda) const
{
    /* Compute the diagonal preconditoned residue : delta = - (D - lambda)^{-1} r */
    int size = r.rows();
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);
    for (int i=0; i < size; i++) {
        delta(i) = r(i) / (lambda - D(i));
    }

    return delta;
}

Eigen::VectorXd DavidsonSolver::_olsen_correction(Eigen::VectorXd &r, Eigen::VectorXd &x, Eigen::VectorXd &D, double lambda) const
{
    /* Compute the olsen correction :

    \delta = (D-\lambda)^{-1} (-r + \epsilon x)

    */

    int size = r.rows();
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);

    delta = DavidsonSolver::_dpr_correction(r,D,lambda);

    double _num = - x.transpose() * delta;
    double _denom = - x.transpose() * DavidsonSolver::_dpr_correction(x,D,lambda);
    double eps = _num / _denom;

    delta += eps * x;

    return delta;
}

Eigen::MatrixXd DavidsonSolver::_QR(Eigen::MatrixXd &A) const
{
    
    int nrows = A.rows();
    int ncols = A.cols();
    ncols = std::min(nrows,ncols);
    
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    return qr.householderQ() * Eigen::MatrixXd::Identity(nrows,ncols);
}

}}