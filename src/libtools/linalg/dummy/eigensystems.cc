/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace votca { namespace tools {

using namespace std;
 
/**
*
* ublas binding for gsl_eigen_symmv
* note that the eigenvalues/eigenvectors are UNSORTED 
* 
*/
bool linalg_eigenvalues_symmetric(const  ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
    throw std::runtime_error("linalg_eigenvalues_symmetric is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL and MKL support");
};


/**
*
* ublas binding for gsl_eigen_symmv
* input matrix type general matrix!
* wrapping gsl_eigen_symmv 
* 
*/
bool linalg_eigenvalues(const ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL or MKL support");
};

/**
*
* ublas binding for gsl_eigen_symmv
* input matrix type general matrix single precision!
* wrapping gsl_eigen_symmv 
* 
*/
bool linalg_eigenvalues( ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL or MKL support");
};

bool linalg_eigenvalues(  ub::vector<float> &E, ub::matrix<float> &V)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL or MKL support");
};


/**
*
* ublas binding for gsl_eigen_symm
* input matrix type general matrix!
* wrapping gsl_eigen_symm leaves input matrix 
* 
*/
bool linalg_eigenvalues( ub::vector<double> &E, ub::matrix<double> &V)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL or MKL support");
};


/*
 * use expert routine to calculate only a subrange of eigenvalues
 */
bool linalg_eigenvalues(ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V , int nmax)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
}

/*
 * use expert routine to calculate only a subrange of eigenvalues single precision
 */
bool linalg_eigenvalues(ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V , int nmax)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
}

bool linalg_eigenvalues_general(const ub::matrix<double> &A,const ub::matrix<double> &B, ub::vector<double> &E, ub::matrix<double> &V)
{
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL and MKL - recompile Votca Tools with GSL or MKL support");
};

}}
