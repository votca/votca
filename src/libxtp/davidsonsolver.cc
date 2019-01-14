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



namespace votca { namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver() : iter_max(1000), tol(1E-6), max_search_space(100) { }
DavidsonSolver::DavidsonSolver(int iter_max) : iter_max(iter_max) , tol(1E-6), max_search_space(100) { }
DavidsonSolver::DavidsonSolver(int iter_max, real_gwbse tol) : iter_max(iter_max), tol(tol), max_search_space(100) { }
DavidsonSolver::DavidsonSolver(int iter_max, real_gwbse tol, int max_search_space) : iter_max(iter_max), tol(tol), max_search_space(max_search_space) { }

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
MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction(OpMat A, VectorXfd u, real_gwbse lambda)
{
    MatrixXfd w;

    // form the projector 
    MatrixXfd P = -u*u.transpose();
    P.diagonal().array() += 1.0;

    // compute the residue
    VectorXfd r = A*u - lambda*u;

    // project the matrix
    // P * (A - lambda*I) * P^T
    MatrixXfd projA = A*P.transpose();
    projA -= lambda*P.transpose();
    projA = P * projA;

    return DavidsonSolver::_solve_linear_system(projA,r);
}

template MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction<MatrixXfd>(MatrixXfd  A, VectorXfd u, real_gwbse lambda);
template MatrixXfd DavidsonSolver::_jacobi_orthogonal_correction<BSE_MF>(BSE_MF A, VectorXfd u, real_gwbse lambda);

template<class OpMat>
void DavidsonSolver::solve(OpMat A, int neigen, int size_initial_guess)
{

    if (this->_debug_)
    {
        cout << endl;
        cout << "===========================" << endl; 
        if(this->jacobi_correction)
        {
            cout << "= Jacobi-Davidson      " <<  endl; 
            cout << "    linsolve : " << this->jacobi_linsolve << endl;
        }
        else
            cout << "= Davidson (DPR)" <<  endl; 
        cout << "===========================" << endl;
        cout << endl;
    }

    double norm;
    int size = A.rows();

    // if argument not provided we default to 0
    // and set to twice the number of eigenvalues
    if (size_initial_guess == 0)
        size_initial_guess = 2*neigen;

    int search_space = size_initial_guess;

    // initialize the guess eigenvector
    VectorXfd Adiag = A.diagonal();    
    MatrixXfd V = DavidsonSolver::_get_initial_eigenvectors(Adiag,size_initial_guess);

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
        cout << "iter\tSearch Space\tNorm" << endl;
    
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

        // diagonalize in the subspace
        // we could replace that with LAPACK ... 
        Eigen::SelfAdjointEigenSolver<MatrixXfd> es(T);
        lambda = es.eigenvalues();
        U = es.eigenvectors();
        
        // Ritz eigenvectors
        q = V.block(0,0,V.rows(),search_space)*U;

        // compute correction vectors
        // and append to V
        for (int j=0; j<size_initial_guess; j++)
        {   

            // jacobi-davidson correction
            if (this->jacobi_correction)
                w = DavidsonSolver::_jacobi_orthogonal_correction<OpMat>(A,q.col(j),lambda(j));
            
            // Davidson DPR
            else  
                w = ( A*q.col(j) - lambda(j)*q.col(j) ) / ( lambda(j) - Adiag(j) );  

            // append the correction vector to the search space
            V.conservativeResize(Eigen::NoChange,V.cols()+1);
            V.col(V.cols()-1) = w;
        }

        // check for convergence
        norm = (lambda.head(neigen) - lambda_old).norm();

        if(_debug_)
            printf("%4d\t%12d\t%4.2e/%.0e\n", iiter,search_space,norm,tol);
        
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
template void DavidsonSolver::solve<BSE_MF>(BSE_MF A, int neigen, int size_initial_guess=0);

}}