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

namespace votca { namespace tools {

using namespace std;

/**
 * ublas binding to GSL  Singular Value Decomposition
 * 
 * A = U S V^T
 * 
 * @param A MxN matrix do decompose. Becomes an MxN orthogonal matrix U
 * @param V NxN orthogonal square matrix
 * @param E NxN diagonal matrix of singular values
 * @return succeeded or not 
 */
bool linalg_singular_value_decomposition( ub::matrix<double> &A, ub::matrix<double> &VT, ub::vector<double> &S )
{
    throw std::runtime_error("linalg_singular_value_decomposition is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
};

}}
