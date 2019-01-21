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

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <chrono>

#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/bse_matrix_free.h>
#include <votca/xtp/davidsonsolver.h>

using boost::format;
using std::flush;

namespace votca { namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver(ctp::Logger &log) : _log(log), iter_max(1000), tol(1E-6), max_search_space(100) { }
// DavidsonSolver::DavidsonSolver(int iter_max) : iter_max(iter_max) , tol(1E-6), max_search_space(100) { }
// DavidsonSolver::DavidsonSolver(int iter_max, real_gwbse tol) : iter_max(iter_max), tol(tol), max_search_space(100) { }
// DavidsonSolver::DavidsonSolver(int iter_max, real_gwbse tol, int max_search_space) : iter_max(iter_max), tol(tol), max_search_space(max_search_space) { }

void DavidsonSolver::set_iter_max(int N) { this->iter_max = N; }
void DavidsonSolver::set_tolerance(real_gwbse eps) { this->tol = eps; }
void DavidsonSolver::set_max_search_space(int N) { this->max_search_space = N;}
void DavidsonSolver::set_jacobi_correction() { this->jacobi_correction = true; }
void DavidsonSolver::set_jacobi_linsolve(int method) {this->jacobi_linsolve = method;}

VectorXfd DavidsonSolver::eigenvalues() {return this->_eigenvalues;}
MatrixXfd DavidsonSolver::eigenvectors() {return this->_eigenvectors;}


ArrayXfd DavidsonSolver::_sort_index(VectorXfd V)
{
    ArrayXfd idx = ArrayXfd::LinSpaced(V.rows(),0,V.rows()-1);
    std::sort(idx.data(),idx.data()+idx.size(),
              [&](int i1, int i2){return V[i1]<V[i2];});
    return idx; 
}

MatrixXfd DavidsonSolver::_get_initial_eigenvectors(VectorXfd d, int size_initial_guess)
{

    MatrixXfd guess = MatrixXfd::Zero(d.size(),size_initial_guess);
    ArrayXfd idx = DavidsonSolver::_sort_index(d);

    for (int j=0; j<size_initial_guess;j++)
        guess(idx(j),j) = 1.0;

    return guess;
}

MatrixXfd DavidsonSolver::_solve_linear_system(MatrixXfd A, VectorXfd r)
{
    MatrixXfd w;

    switch(this->jacobi_linsolve)
    {        
        //use cg approximate solver     
        case 0: 
        {
            Eigen::ConjugateGradient<MatrixXfd, Eigen::Lower|Eigen::Upper> cg;
            cg.compute(A);
            w = cg.solve(r);
        }   

        //use GMRES approximate solver
        case 1:
        {
            Eigen::GMRES<MatrixXfd, Eigen::IdentityPreconditioner> gmres;
            gmres.compute(A);
            w = gmres.solve(r);
        }

        //use llt exact solver
        case 2: w = A.llt().solve(r);
    }

    return w;
}



template <class OpMat>
MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction(OpMat A, VectorXfd r, VectorXfd u, real_gwbse lambda)
{
    MatrixXfd w;

    // form the projector 
    MatrixXfd P = -u*u.transpose();
    P.diagonal().array() += 1.0;

    // compute the residue
    //VectorXfd r = A*u - lambda*u;

    // project the matrix
    // P * (A - lambda*I) * P^T
    MatrixXfd projA = A*P.transpose();
    projA -= lambda*P.transpose();
    projA = P * projA;

    return DavidsonSolver::_solve_linear_system(projA,r);
}

template MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction<MatrixXfd>(MatrixXfd  A, VectorXfd r, VectorXfd u, real_gwbse lambda);
//template MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction<BSE_MF>(BSE_MF A, VectorXfd u, real_gwbse lambda);

template<class OpMat>
void DavidsonSolver::solve(OpMat A, int neigen, int size_initial_guess)
{

    if (this->_debug_)
    {
        
        if(this->jacobi_correction)
            CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Jacobi-Davidson - linsolve : " << this->jacobi_linsolve << flush;
        else
            CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Davidson (DPR)" << flush;
    }

    double norm;
    int size = A.rows();

    // asking too many eigenvalues for the system size
    if (neigen>size/4)
    {
        CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Warning neigen (" << neigen << ") larger than system size (" << size << ")" << flush;
        neigen = size/4;
        CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Computing only " << neigen << " eigenvalues" << flush;                            
    }

    //. search space exeeding the system size
    if (max_search_space > size)
    {
        CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Warning Max search space (" << max_search_space << ") larger than system size (" << size << ")" << flush;   
        max_search_space = size;
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Max search space set to " << size << flush;                            
    }
    //std::cout << std::endl << "matrix A : " << std::endl << A.block(0,0,5,5) << std::endl;

    // if argument not provided we default to 0
    // and set to twice the number of eigenvalues
    if (size_initial_guess == 0)
        size_initial_guess = 2*neigen;

    int search_space = size_initial_guess;

    // initialize the guess eigenvector
    VectorXfd Adiag = A.diagonal();    
    MatrixXfd V = DavidsonSolver::_get_initial_eigenvectors(Adiag,size_initial_guess);
    //std::cout << std::endl << "matrix V : " << std::endl << V << std::endl;


    // sort the eigenvalues
    std::sort(Adiag.data(),Adiag.data()+Adiag.size());

    // thin matrix for QR
    MatrixXfd thinQ;

    // eigenvalues hodlers
    VectorXfd lambda;
    VectorXfd lambda_old = VectorXfd::Ones(neigen,1);

    // temp varialbes 
    MatrixXfd T, U, w, q;

    // chrono !
    chrono::time_point<chrono::system_clock> start, end, instart, instop;
    chrono::duration<double> elapsed_time;

    if (_debug_)
      CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " iter\tSearch Space\tNorm" << flush;
        
    
    for (int iiter = 0; iiter < iter_max; iiter ++ )
    {
        
        // orthogonalise the vectors
        // use the HouseholderQR algorithm of Eigen
        if (iiter>0)
        {
            thinQ = MatrixXfd::Identity(V.rows(),V.cols());
            Eigen::HouseholderQR<MatrixXfd> qr(V);
            V = qr.householderQ() * thinQ;
        }

        // project the matrix on the trial subspace
        T = A * V;
        T = V.transpose()*T;
        //std::cout << std::endl << " T matrix : " << std::endl << T << std::endl;

        // diagonalize in the subspace
        // we could replace that with LAPACK ... 
        Eigen::SelfAdjointEigenSolver<MatrixXfd> es(T);
        lambda = es.eigenvalues();
        U = es.eigenvectors();
        
        // Ritz eigenvectors
        q = V.block(0,0,V.rows(),search_space)*U;

        // compute correction vectors
        // and append to V
        norm = 0.0;
        for (int j=0; j<size_initial_guess; j++)
        {   

            // residue vector
            w = A*q.col(j) - lambda(j)*q.col(j);
            norm += w.norm();

            // jacobi-davidson correction
            if (this->jacobi_correction)
                w = DavidsonSolver::_jacobi_orthogonal_correction<OpMat>(A,w,q.col(j),lambda(j));
            
            // Davidson DPR
            else  
                w = w / ( lambda(j) - Adiag(j) );  

            // append the correction vector to the search space
            V.conservativeResize(Eigen::NoChange,V.cols()+1);
            V.col(V.cols()-1) = w;
        }

        // check for convergence
        //norm = (lambda.head(neigen) - lambda_old).norm();
        norm /= size_initial_guess;

        if(_debug_)
            CTP_LOG(ctp::logDEBUG, _log)  << ctp::TimeStamp() 
               << format(" %1$4d \t %2$12d \t %3$4.2e/%4$.0e") % iiter % search_space % norm % tol << flush; 
            // CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() <<
            //     iiter << "   " << search_space  << "   " << norm << "   " << tol << flush; 
        
        
        // break if converged, update otherwise
        if (norm < tol)
            break;
        else
        {
            lambda_old = lambda.head(neigen);
            search_space += size_initial_guess;
        }

        // restart
        if (search_space > max_search_space)
        {
            V = q.block(0,0,V.rows(),size_initial_guess);
            search_space = size_initial_guess;
        }
    }

    // store the eigenvalues/eigenvectors
    this->_eigenvalues = lambda.head(neigen);
    this->_eigenvectors = U.block(0,0,U.rows(),neigen);

}


template void DavidsonSolver::solve<MatrixXfd>(MatrixXfd A, int neigen, int size_initial_guess=0);
//template void DavidsonSolver::solve<BSE_MF>(BSE_MF &A, int neigen, int size_initial_guess=0);

}}