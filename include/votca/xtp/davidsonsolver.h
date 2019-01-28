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

#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/qmstate.h>

namespace votca { namespace xtp {
 
    /**
    * \brief Use Davidson algorithm to solve A*V=E*V

    **/

    class DavidsonSolver
    {

        public:

            DavidsonSolver(ctp::Logger &log);
            // DavidsonSolver(ctp::Logger &log,int itermax);
            // DavidsonSolver(int itermax, real_gwbse tol);
            // DavidsonSolver(int itermax, real_gwbse tol, int max_search_size);

            void set_iter_max(int set_iter_max);
            void set_tolerance(real_gwbse tol);
            void set_max_search_space(int size);
            void set_jacobi_correction();
            void set_jacobi_linsolve(int method);

            Eigen::VectorXd eigenvalues();
            Eigen::MatrixXd eigenvectors();

            template <typename MatrixReplacement>
            void solve(MatrixReplacement &A, int neigen, int size_initial_guess = 0);

        private :

            ctp::Logger &_log;
            int iter_max;
            double tol;
            int max_search_space;
            bool _debug_ = true;

            bool jacobi_correction=false;
            int jacobi_linsolve = 0;

            Eigen::VectorXd _eigenvalues;
            Eigen::MatrixXd _eigenvectors; 

            Eigen::ArrayXd _sort_index(Eigen::VectorXd &V);
            Eigen::MatrixXd _get_initial_eigenvectors(Eigen::VectorXd &D, int size);
            Eigen::MatrixXd _solve_linear_system(Eigen::MatrixXd &A, Eigen::VectorXd &b); 

            template <typename MatrixReplacement>
            Eigen::MatrixXd _jacobi_orthogonal_correction(MatrixReplacement &A, Eigen::VectorXd &r, Eigen::VectorXd &u, real_gwbse lambda);

    };

   
   
}}



#endif	// __VOTCA_TOOLS_DAVIDSON_SOLVER_H

