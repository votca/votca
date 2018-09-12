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
    
    int size() const{return dimension;}

    double TraceofProd(const Symmetric_Matrix& a) const;

    void AddtoEigenMatrix(Eigen::MatrixXd& full, double factor = 1.0) const;
    
    Eigen::MatrixXd FullMatrix()const;

    double &operator()(const size_t i,const size_t j) {
        return data[Index(i,j)];
    };
    
    
    const double &operator()(const size_t i,const size_t j) const{
        return data[Index(i,j)];
    };
    
    
    friend std::ostream &operator<<(std::ostream &out, const Symmetric_Matrix& a);
private:
    
    size_t Index(const size_t i,const size_t j)const;
       
    std::vector<double> data;
    size_t dimension;
    
  
    
};




 
}}

#endif	/* AOBASIS_H */

