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

            void set_iter_max(int N) { this->iter_max = N; }
            void set_tolerance(double eps) { this->tol = eps; }
            void set_max_search_space(int N) { this->max_search_space = N;}
            void set_jacobi_correction() { this->jacobi_correction = true; }
            void set_jacobi_linsolve(std::string method) {this->jacobi_linsolve = method;}

            Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
            Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }

            template <typename MatrixReplacement>
            void solve(MatrixReplacement &A, int neigen, int size_initial_guess = 0);

        private :

            ctp::Logger &_log;
            int iter_max = 1000;
            double tol = 1E-6;
            int max_search_space = 100;
        
            bool jacobi_correction = false;
            std::string jacobi_linsolve = "CG";

            Eigen::VectorXd _eigenvalues;
            Eigen::MatrixXd _eigenvectors; 

            Eigen::ArrayXd _sort_index(Eigen::VectorXd &V) const;
            Eigen::MatrixXd _get_initial_eigenvectors(Eigen::VectorXd &D, int size) const;
            Eigen::MatrixXd _solve_linear_system(Eigen::MatrixXd &A, Eigen::VectorXd &b) const; 

            template <typename MatrixReplacement>
            Eigen::MatrixXd _jacobi_orthogonal_correction(MatrixReplacement &A, Eigen::VectorXd &r, Eigen::VectorXd &u, double lambda) const;

    };

   
   
}}



#endif	// __VOTCA_TOOLS_DAVIDSON_SOLVER_H

