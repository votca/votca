/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include "votca/xtp/symmetric_matrix.h"
#include <iostream>



namespace votca { namespace xtp {
 
  
  Symmetric_Matrix::Symmetric_Matrix(const Eigen::MatrixXd& full):dimension(full.rows()) {
        assert(full.rows() == full.cols() && "Input matrix not quadratic");
       data.resize((dimension + 1) * dimension / 2);
        for (int i = 0; i < full.rows(); ++i) {
            for (int j = 0; j <= i; ++j) {
                this->operator()(i, j) = full(i, j);
                //std::cout<<i<<" "<<j<<" "<<full(i, j)<<" "<<this->operator()(i, j)<<std::endl;
            }
        }
    }
  
   std::ostream &operator<<(std::ostream &out, const Symmetric_Matrix& a) {

    out << "[" << a.dimension << "," << a.dimension << "]\n";
    for (unsigned i = 0; i < a.dimension; ++i) {
      for (unsigned j = 0; j <= i; ++j) {
        out<<a(i,j);
        if(i==j){
          out<<"\n";
        }else{
          out<<" ";
        }
            
     }
    }
            return out;
        }

}}
