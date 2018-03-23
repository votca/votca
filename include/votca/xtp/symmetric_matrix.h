/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef __XTP_SYMMETRIC_MATRIX__H
#define	__XTP_SYMMETRIC_MATRIX__H



#include <votca/xtp/basisset.h>





namespace votca { namespace xtp {

/*
 * A symmetric matrix implementation for doubles, acces upper diagonal matrix
 */
class Symmetric_Matrix
{   
public:

Symmetric_Matrix(size_t dim) {
            dimension = dim;
            data.resize((dim + 1) * dim / 2);
        }

        int size() {
            return size();
        }

    Symmetric_Matrix(Eigen::MatrixXd full) {
        assert(full.rows() == full.cols() && "Input matrix not quadratic");
        dimension = full.rows();
        for (int i = 0; i < full.rows(); ++i) {
            for (int j = 0; j <= i; ++j) {
                this->operator(i, j) = full(i, j);
            }
        }
    }

    Symmetric_Matrix operator+(const Symmetric_Matrix& a, const Symmetric_Matrix& b) {
        Symmetric_Matrix c(a.dimension);
        for (size_t i = 0; i < dimension; ++i) {
            c.data[i] = a.data[i] + b.data[i];
        }
        return c;
    }

    void AddtoEigenMatrix(Eigen::MatrixXd& full, double factor = 1.0) {
        for (int i = 0; i < full.rows(); ++i) {
            for (int j = 0; j < full.cols(); ++j) {
                full(i, j) = factor * this->operator(i, j);
            }
        }
        return;
    }

    double &operator()(size_t i, size_t j) {
        size_t index;

        if (i < j) {
            index = i * (i + 1) / 2 + j;
        } else {
            index = j * (j + 1) / 2 + i;
        }
        return data[index];
    };
    
private:
       
    std::vector<double> data;
    size_t dimension;
    
  
    
};


 
}}

#endif	/* AOBASIS_H */

