/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <fstream>
#include <math.h>
#include <votca/tools/property.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>

namespace votca {
namespace xtp {

void PadeApprox::clear() {

  _value.clear();

  _grid.clear();

  _coeff.clear();
  
  _rejected_points=0;
}

void PadeApprox::addPoint(std::complex<double> frequency, Eigen::Matrix3d value) {

  this->_grid.push_back(frequency);
  this->_value.push_back(value);
  this->_coeff.push_back(RecursivePolynom(_grid.size() - 1, _grid.size()));
  if (_coeff.at(_coeff.size() - 1).norm() !=
      _coeff.at(_coeff.size() - 1).norm()) {
    std::cout
        << "reject point, unvalid coeficient at w=" << frequency
        << std::endl;  //<<_coeff.at(_coeff.size()-1)<<std::endl<<std::endl;
    _coeff.pop_back();
    _grid.pop_back();
    _value.pop_back();
    _rejected_points++;
  }
}

Eigen::Matrix3d PadeApprox::RecursivePolynom(int indx, int degree) {
    
  if (degree == 1) {
    //std::cout<<"g_"<<degree<<"("<<_grid.at(indx)<<")="<<_value.at(indx)(0,0)<<std::endl;
    return _value.at(indx);
  } else if(this->_temp_container_g[indx][degree-1](0,0)!=0){
      return this->_temp_container_g[indx][degree-1];
  }else{
    Eigen::Matrix3d temp = RecursivePolynom(indx, degree - 1);
    Eigen::Matrix3d u;
    if(indx==degree-1){
        //std::cout<<"Used old coefficient for g with "<<degree-1<<" and index "<<indx-1<<std::endl;
        //std::cout<<"Coeff = "<<_coeff[indx-1](0,0)<<std::endl;
        u = _coeff[indx-1] - temp;
        //u = RecursivePolynom(indx-1, degree-1) - temp;
    } else{
        u = RecursivePolynom(degree-2, degree-1) - temp;
    }
    
    //std::cout<<"u="<<u(0,0)<<std::endl;
    
    Eigen::Matrix3d l = (temp * (_grid.at(indx) - _grid.at(degree - 2))).real();
    
    //std::cout<<"l="<<l(0,0)<<std::endl;
    //std::cout<<"u= "<<u<<std::endl;
    //std::cout<<"l inverse check "<<std::endl<<l*l.inverse()<<std::endl<<std::endl;
    //std::cout<<"det(l)"<<l.determinant()<<std::endl;
    Eigen::Matrix3d inverse;
    double determinant;
    bool invertible;
    l.computeInverseAndDetWithCheck(inverse,determinant,invertible);
    if(!invertible){
        std::cout<<"Matrix not invertible "<<determinant<<std::endl<<"l="<<std::endl<<l<<std::endl<<"inverse"<<std::endl<<inverse<<std::endl<<"Product"<<std::endl<<l*inverse<<std::endl;
    }
    std::cout<<"l inverse check "<<std::endl<<l*inverse<<std::endl<<std::endl;
    Eigen::Matrix3d result=u * (inverse);
    _temp_container_g[indx][degree-1]=result;
    return result;
  }
}

Eigen::Matrix3d PadeApprox::RecursiveA(std::complex<double> frequency, int index) {

    if(_temp_container_A.size()>index){
        return _temp_container_A[index];
    }
    else {
      Eigen::Matrix3d A= ((frequency - _grid.at(index - 2)) *
                                                    _coeff.at(index - 1) *
                                                    RecursiveA(frequency, index - 2) + RecursiveA(frequency, index - 1)).real();
    _temp_container_A.push_back(A);
    return A;
  }
}

Eigen::Matrix3d PadeApprox::RecursiveB(std::complex<double> frequency, int index) {

  if(_temp_container_B.size()>index){
        return _temp_container_B[index];
  }
  else {
    Eigen::Matrix3d B= ((frequency - _grid.at(index - 2)) *
                                                    _coeff.at(index - 1) *
                                                    RecursiveB(frequency, index - 2) + RecursiveB(frequency, index - 1)).real();
    _temp_container_B.push_back(B);
    return B;
  }
}

void PadeApprox::initialize(int num_points) { 
    this->_temp_container_g.resize(num_points);
    for (int i = 0; i < num_points; i++){
        _temp_container_g[i].resize(num_points);
    }
}

Eigen::Matrix3d PadeApprox::evaluatePoint(std::complex<double> frequency) {
    _temp_container_A.clear();
    _temp_container_B.clear();
    _temp_container_A.push_back(Eigen::Matrix3d::Zero());
    _temp_container_A.push_back(_coeff.at(0));
    _temp_container_B.push_back(Eigen::Matrix3d::Identity());
    _temp_container_B.push_back(Eigen::Matrix3d::Identity());
  return RecursiveB(frequency, _grid.size()).inverse() *
         RecursiveA(frequency, _grid.size());
}

}  // namespace xtp
}  // namespace votca