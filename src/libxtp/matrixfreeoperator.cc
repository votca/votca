
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


    MatrixFreeOperator::MatrixFreeOperator(){}

    // virtual here : get a row of the operator
    RowVectorXfd MatrixFreeOperator::row(int index) const
    {
    	
        std::cout << "MatrixFreeOperator.row() not defined in class" << std::endl;
        exit(0);
    }

    // virtual here : get a col of the operator
    VectorXfd MatrixFreeOperator::col(int index) const
    {

        std::cout << "MatrixFreeOperator.col() not defined in class" << std::endl;
        exit(0);
    }

    VectorXfd MatrixFreeOperator::diagonal()
    {
        VectorXfd D = VectorXfd::Zero(OpSizeVal,1);
        RowVectorXfd row_data;
        for(int i=0; i<OpSizeVal;i++)
        {
            row_data = this->row(i);
            D(i,0) = row_data(0,i);
        }
        return D;
    }

    // get the full matrix if we have to
    MatrixXfd MatrixFreeOperator::get_full_mat()
    {
    	MatrixXfd matrix = MatrixXfd::Zero(OpSizeVal,OpSizeVal);
        for(int i=0; i<OpSizeVal; i++)
            matrix.row(i) = this->row(i);
        return matrix; 
    }


    // get the size
    int MatrixFreeOperator::OpSize()
    {
    	return OpSizeVal;
    }

}}