/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/tools/votca_config.h>

#ifndef NOMKL
#include "mkl.h"
#include "mkl_lapacke.h"
#endif


namespace votca { namespace tools {

using namespace std;


bool linalg_singular_value_decomposition( ub::matrix<double> &A, ub::matrix<double> &V, ub::vector<double> &S ){
    
#ifdef NOMKL
    throw std::runtime_error("linalg_singular_value_decomposition is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
#else
        // matrix inversion using MKL
   
    throw std::runtime_error("linalg_singular_value_decomposition is not yet implemented with MKL support");
    return false;
#endif   
}

}}
