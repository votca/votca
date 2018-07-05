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

#ifndef __XTP_SYMMETRIC_MATRIX__H
#define	__XTP_SYMMETRIC_MATRIX__H



#include <votca/xtp/eigen.h>
#include <iostream>
#include <vector>





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

        
    Symmetric_Matrix(const Eigen::MatrixXd& full);
    
    int size() {
            return dimension;
        }

    double TraceofProd(const Symmetric_Matrix& a) const{
        assert(data.size()==a.data.size()&&"Matrices do not have same size");
        double result=0.0;
        
        for (size_t i=0;i<dimension;++i){
            result+=this->operator ()(i,i)*a.operator ()(i,i);
        }
        
        for (size_t i=0;i<dimension;++i){
            for (size_t j=0;j<i;++j){
                result+=2*this->operator ()(i,j)*a.operator ()(i,j);
        }
        }
        return result;
    }

    void AddtoEigenMatrix(Eigen::MatrixXd& full, double factor = 1.0) const{
        for (int i = 0; i < full.rows(); ++i) {
            for (int j = 0; j < full.cols(); ++j) {
                full(i, j) += factor * this->operator ()(i,j);
            }
        }
        return;
    }
    
    Eigen::MatrixXd FullMatrix(){
        Eigen::MatrixXd result=Eigen::MatrixXd(dimension,dimension);
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                result(i, j) =this->operator ()(i,j);
            }
        }
        return result;
    }

    double &operator()(const size_t i,const size_t j) {
        size_t index;

        if (i >= j) {
            index = (i * (i + 1)) / 2 + j;
        } else {
            index = (j * (j + 1)) / 2 + i;
        }

        return data[index];
    };
    
    
    const double &operator()(const size_t i,const size_t j) const{
        size_t index;

        if (i >= j) {
            index = (i * (i + 1)) / 2 + j;
        } else {
            index = (j * (j + 1)) / 2 + i;
        }
        return data[index];
    };
    
    
    friend std::ostream &operator<<(std::ostream &out, const Symmetric_Matrix& a);
private:
       
    std::vector<double> data;
    size_t dimension;
    
  
    
};




 
}}

#endif	/* AOBASIS_H */

