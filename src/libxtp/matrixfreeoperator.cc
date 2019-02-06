
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
#include <iostream>
#include <votca/xtp/matrixfreeoperator.h>


namespace votca { namespace xtp {


    MatrixFreeOperator::MatrixFreeOperator() {}

    // virtual here : get a col of the operator
    Eigen::VectorXd MatrixFreeOperator::col(int index) const 
    {
         throw std::runtime_error("MatrixFreeOperator.col() not defined in class");
    }

    Eigen::VectorXd MatrixFreeOperator::diagonal() const
    {
        Eigen::VectorXd D = Eigen::VectorXd::Zero(_size,1);
        Eigen::VectorXd col_data;
        for(int i=0; i<_size;i++)
        {
            col_data = this->col(i);
            D(i,0) = col_data(i,0);
        }
        return D;
    }

    // get the full matrix if we have to
    Eigen::MatrixXd MatrixFreeOperator::get_full_matrix() const
    {
        std::cout << " Matrix Size : " << _size << "x" << _size ;
    	Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(_size,_size);
        for(int i=0; i<_size; i++)
            matrix.col(i) = this->col(i);
        std::cout << " done " << std::endl;
        return matrix; 
    }


    // get the size
    int MatrixFreeOperator::size() const
    {
    	return this->_size;
    }

    // set the size
    void MatrixFreeOperator::set_size(int size)
    {
        std::cout << "Set Matrix Size : " << size << "x" << size << std::endl;
        this->_size = size;
    }

}}