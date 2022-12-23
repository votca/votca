/* 
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Standard library includes
#include <iostream>

// VOTCA includes
#include <votca/tools/eigen.h>

// Local VOTCA includes
#include "votca/csg/ml.h"

namespace votca { namespace csg {

using namespace votca::tools;
using namespace std;
  
    
ML::ML()
{
    //initialize parameters with some value
    _bead_number=0;
    _struct_number=0;
    _container_number=0;
    _descriptor_size=0;
    //default is _threebody=0
    _threebody=0;
    for (int i=0; i<3; ++i){
        _sigma(i)=0;
    }

    //initialize as zero vectors/matrix with length zero
    _descriptor=Eigen::VectorXd::Zero(0);
    _descriptor_number=Eigen::VectorXi::Zero(0);
    _alpha=Eigen::VectorXd::Zero(0);
    _L=Eigen::MatrixXd::Zero(0,0);
}

ML::~ML()
{
    //TODO 
}

ML::ML(int bead_number, int struct_number, Eigen::Vector3d sigma, bool threebody)//, double lambda)
{
    //set parameters
    _bead_number=bead_number;
    _struct_number=struct_number;
    _container_number=0;
    _threebody=threebody;
    //descriptor size depends on threebody or not
    if (_threebody){
        _descriptor_size=_struct_number*9;
    } else {
        _descriptor_size=_struct_number*3;
    }
    _descriptor_size=_struct_number*9;
    //initialize with values of sigma
    for (int i=0; i<3; ++i){
        _sigma(i)=sigma(i);
    }
    //_lambda=lambda;

    //initialize corresponding vectors and matrices
    _descriptor=Eigen::VectorXd::Zero(_descriptor_size);
    _descriptor_number=Eigen::VectorXi::Zero(_struct_number);
    _alpha=Eigen::VectorXd::Zero(_bead_number*3);
    _L=Eigen::MatrixXd::Zero(_bead_number*3,_struct_number*3);
}

void ML::Resize(int bead_number, int struct_number){

    //temporarily store the old parameter values
    int bead_number_old = _bead_number;
    int struct_number_old = _struct_number;
    int descriptor_size_old = _descriptor_size;

    //set parameters to the new values
    _bead_number=bead_number;
    _struct_number=struct_number;
    //descriptor size depends on threebody or not
    if (_threebody){
        _descriptor_size=_struct_number*9;
    } else {
        _descriptor_size=_struct_number*3;
    }

    //resize corresponding vectors and matrices, preserving the existing entries
    //the new entries are not initialized, and have to be set to zero
    _descriptor.conservativeResize(_descriptor_size);
    for (int i=descriptor_size_old; i<_descriptor_size; ++i){
        _descriptor(i)=0;
    }
    _descriptor_number.conservativeResize(_struct_number);
    for (int i=struct_number_old; i<_struct_number; ++i){
        _descriptor_number(i)=0;
    }
    _alpha.conservativeResize(_bead_number*3);
    for (int i=bead_number_old*3; i<_bead_number*3; ++i){
        _alpha(i)=0;
    }
    _L.conservativeResize(_bead_number*3,_struct_number*3);
    for (int i=0; i<_bead_number*3; ++i){
        for (int j=0; j<_struct_number*3; ++j){
            //new entries correspond to either i>=bead_number_old or j>=struct_number_old
            if ( (i >= bead_number_old*3) || (j >=struct_number_old*3) ){
                _L(i,j)=0.0;
            }
        }
    }
}

const Eigen::MatrixXd ML::Evaluate_L_K_L_T()
{
    //matrix to store Kernel matrix
    Eigen::MatrixXd K;
    K = Eigen::MatrixXd::Zero(_struct_number*3,_struct_number*3);
    
    //descriptor size depends on threebody or not
    if (_threebody){
        //evaluate Kernel matrix
        EvaluateKernel_Threebody(_descriptor,_descriptor_number,_descriptor,_descriptor_number,K);
    } else {
        //evaluate Kernel matrix
        EvaluateKernel(_descriptor,_descriptor,K);
    }

    //matrix to store product L*K*L_T
    Eigen::MatrixXd L_K_L_T;
    L_K_L_T = Eigen::MatrixXd::Zero(_bead_number*3,_bead_number*3);

    L_K_L_T = _L*K*_L.transpose();

    return L_K_L_T;
}

void ML::PredictStruct(Eigen::VectorXd &descriptor_pred, Eigen::VectorXd &pred, Eigen::MatrixXd &L_K_pred)
{
    //Kernel matrix between descriptor_pred and _descriptor of training data
    Eigen::MatrixXd K_pred;
    K_pred = Eigen::MatrixXd::Zero(_struct_number*3,descriptor_pred.size());
    //Evaluate Kernel matrix
    EvaluateKernel(_descriptor,descriptor_pred,K_pred);
    //matrix to store product L*K_pred
    L_K_pred = _L*K_pred;
    //now evaluate prediction (is equal to _alpha * _L*K_pred)
    pred = L_K_pred.transpose()*_alpha;
}

void ML::PredictStruct_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::VectorXd &pred, Eigen::MatrixXd &L_K_pred)
{
    //Kernel matrix between descriptor_pred and _descriptor of training data
    Eigen::MatrixXd K_pred;
    K_pred = Eigen::MatrixXd::Zero(_struct_number*3,descriptor_pred.size()/3);
    //Evaluate Kernel matrix
    EvaluateKernel_Threebody(_descriptor,_descriptor_number,descriptor_pred,descriptor_pred_number,K_pred);
    //matrix to store product L*K_pred
    L_K_pred = _L*K_pred;
    //now evaluate prediction (is equal to _alpha * _L*K_pred)
    pred = L_K_pred.transpose()*_alpha;
}

void ML::PredictStruct(Eigen::VectorXd &descriptor_pred, Eigen::VectorXd &pred)
{
    //Kernel matrix between descriptor_pred and _descriptor of training data
    Eigen::MatrixXd K_pred;
    K_pred = Eigen::MatrixXd::Zero(_struct_number*3,descriptor_pred.size());
    //Evaluate Kernel matrix
    EvaluateKernel(_descriptor,descriptor_pred,K_pred);
    Eigen::MatrixXd L_K_pred;
    L_K_pred = Eigen::MatrixXd::Zero(_bead_number*3, descriptor_pred.size());
    //matrix to store product L*K_pred
    L_K_pred = _L*K_pred;
    //now evaluate prediction (is equal to _alpha * _L*K_pred)
    pred = L_K_pred.transpose()*_alpha;
}

void ML::PredictStruct_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::VectorXd &pred)
{
    //Kernel matrix between descriptor_pred and _descriptor of training data
    Eigen::MatrixXd K_pred;
    K_pred = Eigen::MatrixXd::Zero(_struct_number*3,descriptor_pred.size()/3);
    //Evaluate Kernel matrix
    EvaluateKernel_Threebody(_descriptor,_descriptor_number,descriptor_pred,descriptor_pred_number,K_pred);
    Eigen::MatrixXd L_K_pred;
    L_K_pred = Eigen::MatrixXd::Zero(_bead_number*3, descriptor_pred.size()/3);
    //matrix to store product L*K_pred
    L_K_pred = _L*K_pred;
    //now evaluate prediction (is equal to _alpha * _L*K_pred)
    pred = L_K_pred.transpose()*_alpha;
}

void ML::PredictBead(Eigen::VectorXd &descriptor_pred, Eigen::MatrixXd &mapping_pred, Eigen::VectorXd &pred_beads)
{
    //to store prediction of the structural units
    Eigen::VectorXd pred;
    pred = Eigen::VectorXd::Zero(descriptor_pred.size());
    //prediction of values for structural units
    Eigen::MatrixXd L_K_pred;
    L_K_pred = Eigen::MatrixXd::Zero(_bead_number*3, descriptor_pred.size());
    PredictStruct(descriptor_pred, pred, L_K_pred);
    //now evaluate prediction of the beads (is equal to L_pred * (K_pred * _alpha)) (Vector storing the ML coefficients)
    pred_beads=mapping_pred*pred;
}

void ML::PredictBead_Threebody(Eigen::VectorXd &descriptor_pred, Eigen::VectorXi &descriptor_pred_number, Eigen::MatrixXd &mapping_pred, Eigen::VectorXd &pred_beads)
{
    //to store prediction of the structural units
    Eigen::VectorXd pred;
    pred = Eigen::VectorXd::Zero(descriptor_pred.size()/3);
    //prediction of values for structural units
    Eigen::MatrixXd L_K_pred;
    L_K_pred = Eigen::MatrixXd::Zero(_bead_number*3, descriptor_pred.size()/3);
    PredictStruct_Threebody(descriptor_pred, descriptor_pred_number, pred, L_K_pred);
    //now evaluate prediction of the beads (is equal to L_pred * (K_pred * _alpha)) (Vector storing the ML coefficients)
    pred_beads=mapping_pred*pred;
}

void ML::EvaluateKernel(Eigen::VectorXd &descriptor1, Eigen::VectorXd &descriptor2, Eigen::MatrixXd &K)
{
    //construct kernel matrix
    for (int i=0; i<descriptor1.size()/3.0; i++){
        for (int j=0; j<descriptor2.size()/3.0; j++){
            Eigen::Vector3d vec1(descriptor1(i*3), descriptor1(i*3+1), descriptor1(i*3+2));
            Eigen::Vector3d vec2(descriptor2(j*3), descriptor2(j*3+1), descriptor2(j*3+2));
            double e;
            e = exp( ( -0.5*( (vec1.norm()-vec2.norm())*(vec1.norm()-vec2.norm()) ) )/ (_sigma(0)*_sigma(0)) )*(1.0/(_sigma(0)*_sigma(0)))*(1-( (vec1.norm()-vec2.norm())*(vec1.norm()-vec2.norm()) )/ (_sigma(0)*_sigma(0)) );
            Eigen::Vector3d v1;
            Eigen::Vector3d v2;
            v1 = vec1.normalized();
            v2 = vec2.normalized();
            Eigen::Matrix3d tK = e * v1 * v2.transpose();
            for (int k=0; k<3; ++k){
                for (int l=0; l<3; ++l){
                    K(i*3+k,j*3+l) = tK(k,l);
                }
            }
        }
    }
}

void ML::EvaluateKernel_Threebody(Eigen::VectorXd &descriptor1, Eigen::VectorXi &descriptor_number1, Eigen::VectorXd &descriptor2, Eigen::VectorXi &descriptor_number2, Eigen::MatrixXd &K)
{
    //construct kernel matrix 
    for (int i=0; i<descriptor1.size()/9.0; i++){
        for (int j=0; j<descriptor2.size()/9.0; j++){
            Eigen::Vector3d vec1_12(descriptor1(i*9+0), descriptor1(i*9+1), descriptor1(i*9+2));
            Eigen::Vector3d vec1_13(descriptor1(i*9+3), descriptor1(i*9+4), descriptor1(i*9+5));
            Eigen::Vector3d vec1_23(descriptor1(i*9+6), descriptor1(i*9+7), descriptor1(i*9+8));
            Eigen::Vector3d vec2_12(descriptor2(j*9+0), descriptor2(j*9+1), descriptor2(j*9+2));
            Eigen::Vector3d vec2_13(descriptor2(j*9+3), descriptor2(j*9+4), descriptor2(j*9+5));
            Eigen::Vector3d vec2_23(descriptor2(j*9+6), descriptor2(j*9+7), descriptor2(j*9+8));

            Eigen::Vector3d var1, var2;
            double e;
            var1.x() = vec1_12.norm();
            var1.y() = vec1_13.norm();
            var1.z() = vec1_23.norm();

            var2.x() = vec2_12.norm();
            var2.y() = vec2_13.norm();
            var2.z() = vec2_23.norm();

            //now e with individual sigma values
            e = exp(-0.5*( ( (var1.x()-var2.x())*(var1.x()-var2.x())/(_sigma(0)*_sigma(0)) )+( (var1.y()-var2.y())*(var1.y()-var2.y())/(_sigma(1)*_sigma(1)) )+( (var1.z()-var2.z())*(var1.z()-var2.z())/(_sigma(2)*_sigma(2)) ) ) );

            double K_der_11 = e*( (1/(_sigma(0)*_sigma(0))) - ( (1/(_sigma(0)*_sigma(0)))*(var1.x()-var2.x()) )*( (1/(_sigma(0)*_sigma(0)))*(var1.x()-var2.x()) ) );
            double K_der_22 = e*( (1/(_sigma(1)*_sigma(1))) - ( (1/(_sigma(1)*_sigma(1)))*(var1.y()-var2.y()) )*( (1/(_sigma(1)*_sigma(1)))*(var1.y()-var2.y()) ) );
            double K_der_33 = e*( (1/(_sigma(2)*_sigma(2))) - ( (1/(_sigma(2)*_sigma(2)))*(var1.z()-var2.z()) )*( (1/(_sigma(2)*_sigma(2)))*(var1.z()-var2.z()) ) );

            Eigen::Vector3d vec_q_1_11,vec_q_1_21,vec_q_1_31,vec_q_1_12,vec_q_1_22,vec_q_1_32,vec_q_1_13,vec_q_1_23,vec_q_1_33;

            vec_q_1_11 = (-(vec1_12)/(vec1_12.norm()));
            vec_q_1_12 = ((vec1_12)/(vec1_12.norm()));
            vec_q_1_13(0) = 0;
            vec_q_1_13(1) = 0;
            vec_q_1_13(2) = 0;
            vec_q_1_21 = (-(vec1_13)/(vec1_13.norm()));
            vec_q_1_22(0) = 0;
            vec_q_1_22(1) = 0;
            vec_q_1_22(2) = 0;
            vec_q_1_23 = ((vec1_13)/(vec1_13.norm()));
            vec_q_1_31(0) = 0;
            vec_q_1_31(1) = 0;
            vec_q_1_31(2) = 0;
            vec_q_1_32 = (-(vec1_23)/(vec1_23.norm()));
            vec_q_1_33 = ((vec1_23)/(vec1_23.norm()));

            Eigen::Vector3d vec_q_2_11,vec_q_2_21,vec_q_2_31,vec_q_2_12,vec_q_2_22,vec_q_2_32,vec_q_2_13,vec_q_2_23,vec_q_2_33;

            vec_q_2_11 = (-(vec2_12)/(vec2_12.norm()));
            vec_q_2_12 = ((vec2_12)/(vec2_12.norm()));
            vec_q_2_13(0) = 0;
            vec_q_2_13(1) = 0;
            vec_q_2_13(2) = 0;
            vec_q_2_21 = (-(vec2_13)/(vec2_13.norm()));
            vec_q_2_22(0) = 0;
            vec_q_2_22(1) = 0;
            vec_q_2_22(2) = 0;
            vec_q_2_23 = ((vec2_13)/(vec2_13.norm()));
            vec_q_2_31(0) = 0;
            vec_q_2_31(1) = 0;
            vec_q_2_31(2) = 0;
            vec_q_2_32 = (-(vec2_23)/(vec2_23.norm()));
            vec_q_2_33 = ((vec2_23)/(vec2_23.norm()));

            Eigen::Matrix3d tK;

            if (descriptor_number1(i) == 1){
                if (descriptor_number2(j) == 1){
                    tK =  K_der_11 * vec_q_1_11 * vec_q_2_11.transpose() + K_der_22 * vec_q_1_21 * vec_q_2_21.transpose() + K_der_33 * vec_q_1_31 * vec_q_2_31.transpose();
                }
                if (descriptor_number2(j) == 2){
                    tK =  K_der_11 * vec_q_1_11 * vec_q_2_12.transpose() + K_der_22 * vec_q_1_21 * vec_q_2_22.transpose() + K_der_33 * vec_q_1_31 * vec_q_2_32.transpose();
                }
                if (descriptor_number2(j) == 3){
                    tK =  K_der_11 * vec_q_1_11 * vec_q_2_13.transpose() + K_der_22 * vec_q_1_21 * vec_q_2_23.transpose() + K_der_33 * vec_q_1_31 * vec_q_2_33.transpose();
                }
            }
            if (descriptor_number1(i) == 2){
                if (descriptor_number2(j) == 1){
                    tK =  K_der_11 * vec_q_1_12 * vec_q_2_11.transpose() + K_der_22 * vec_q_1_22 * vec_q_2_21.transpose() + K_der_33 * vec_q_1_32 * vec_q_2_31.transpose();
                }
                if (descriptor_number2(j) == 2){
                    tK =  K_der_11 * vec_q_1_12 * vec_q_2_12.transpose() + K_der_22 * vec_q_1_22 * vec_q_2_22.transpose() + K_der_33 * vec_q_1_32 * vec_q_2_32.transpose();
                }
                if (descriptor_number2(j) == 3){
                    tK =  K_der_11 * vec_q_1_12 * vec_q_2_13.transpose() + K_der_22 * vec_q_1_22 * vec_q_2_23.transpose() + K_der_33 * vec_q_1_32 * vec_q_2_33.transpose();
                }
            }
            if (descriptor_number1(i) == 3){
                if (descriptor_number2(j) == 1){
                    tK =  K_der_11 * vec_q_1_13 * vec_q_2_11.transpose() + K_der_22 * vec_q_1_23 * vec_q_2_21.transpose() + K_der_33 * vec_q_1_33 * vec_q_2_31.transpose();
                }
                if (descriptor_number2(j) == 2){
                    tK =  K_der_11 * vec_q_1_13 * vec_q_2_12.transpose() + K_der_22 * vec_q_1_23 * vec_q_2_22.transpose() + K_der_33 * vec_q_1_33 * vec_q_2_32.transpose();
                }
                if (descriptor_number2(j) == 3){
                    tK =  K_der_11 * vec_q_1_13 * vec_q_2_13.transpose() + K_der_22 * vec_q_1_23 * vec_q_2_23.transpose() + K_der_33 * vec_q_1_33 * vec_q_2_33.transpose();
                }
            }

            for (int k=0; k<3; ++k){
                for (int l=0; l<3; ++l){
                    K(i*3+k,j*3+l) = tK(k,l);
                }
            }
        }
    }
}

}}


