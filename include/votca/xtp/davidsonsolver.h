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

#ifndef __VOTCA_TOOLS_DAVIDSON_SOLVER_H
#define	__VOTCA_TOOLS_DAVIDSON_SOLVER_H
#include <votca/xtp/eigen.h>

namespace votca { namespace xtp {
 
    /**
    * \brief Use Davidson algorithm to solve A*V=E*V

    **/

    class DavidsonSolver
    {

        public:

            DavidsonSolver();
            DavidsonSolver(int itermax);
            DavidsonSolver(int itermax, real_gwbse tol);
            DavidsonSolver(int itermax, real_gwbse tol, int max_search_size);

            void set_iter_max(int set_iter_max);
            void set_tolerance(real_gwbse tol);
            void set_max_search_space(int size);
            void set_jacobi_correction();
            void set_jacobi_linsolve(int method);

            VectorXfd eigenvalues();
            MatrixXfd eigenvectors();

            template <typename OpMat>
            void solve(OpMat A, int neigen, int size_initial_guess = 0);

        private :

            int iter_max;
            double tol;
            int max_search_space;
            bool _debug_ = true;

            bool jacobi_correction=false;
            int jacobi_linsolve = 0;

            VectorXfd _eigenvalues;
            MatrixXfd _eigenvectors; 

            ArrayXfd _sort_index(VectorXfd V);
            MatrixXfd _get_initial_eigenvectors(VectorXfd D, int size);
            MatrixXfd _solve_linear_system(MatrixXfd, VectorXfd b); 

            template <typename OpMat>
            MatrixXfd _jacobi_orthogonal_correction(OpMat A, VectorXfd u, real_gwbse lambda);

    };

   
   
}}



#endif	// __VOTCA_TOOLS_DAVIDSON_SOLVER_H

