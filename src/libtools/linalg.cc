/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/tools/linalg.h>

#include <iostream>
#include <sstream>
namespace votca { namespace tools {
  
void linalg_constrained_qrsolve(Eigen::VectorXd &x, Eigen::MatrixXd &A, const Eigen::VectorXd &b,const Eigen::MatrixXd &constr)
{
    // check matrix for zero column
    int nonzero_found = 0;
    for(size_t j=0; j<A.cols(); j++) {
        nonzero_found = A.col(j).isApproxToConstant(0.0,1e-9);
        if(nonzero_found==0) {
            throw std::runtime_error("constrained_qrsolve_zero_column_in_matrix");
        }
    }

    const int N = b.size();
    const int ngrid = x.size()/2;
    const int ysize = constr.rows(); //number of constraints is number of rows of constr
    
    
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(constr.transpose());
    Eigen::VectorXd tau_qr=QR.hCoeffs();

    // temporary variables
    Eigen::MatrixXd Q=Eigen::MatrixXd::Identity(2*ngrid, 2*ngrid);  
    // Q matrix: QR decomposition of trans(B) (trans(B) has 2*ngrid rows and ysize columns )
    // Q matrix is orthogonal 2*ngrid x 2*ngrid matrix
    Eigen::VectorXd v=Eigen::VectorXd::Zero(2*ngrid);
    

    //construct the matrix Q of trans(B)=Q * R
    //k=1,..,ysize (number of rows of B, meaning number of constraints)
    for (int k = ysize; k > 0 ; k--) {

        for (int icout = 0; icout < k - 1; icout++) {
             v(icout) = 0;
        }
        v(k - 1) = 1.0;

        for (int icout = k; icout < 2*ngrid; icout++) {
             v(icout) = constr(k - 1,icout );
        }

        Eigen::MatrixXd Q_k = Eigen::MatrixXd::Identity(2*ngrid, 2*ngrid) - tau_qr(k - 1 ) *  (v*v.transpose()) ;
        Q = Q* Q_k;
    }

    //now Q is the one of trans(B)=Q * R
    // Calculate A * Q and store the result in A
    A = A*Q.transpose();


    // A = [A1 A2], so A2 is just a block of A
    // [A1 A2] has N rows. A1 has ysize columns 
    //A2 has 2*ngrid-ysize columns 
    Eigen::MatrixXd A2 = A.block(0,ysize,N,2*ngrid-ysize);
    //now perform QR-decomposition of A2 to solve the least-squares problem A2 * z = b
    //A2 has N rows and (2*ngrid-ysize) columns -> 
    Eigen::VectorXd z=A2.householderQr().solve(b);

    
    // Next two steps assemble vector from y (which is zero-vector) and z
    // (which we just got by gsl_linalg_QR_lssolve)

    Eigen::VectorXd result=Eigen::VectorXd::Zero(2*ngrid);
    

    //vector z has (2*ngrid-ysize) components
    for (int i = ysize; i < 2 * ngrid; i++ ) {
           result[i] = z(i - ysize);
    }
    
    // To get the final answer this vector should be multiplied by matrix Q
    // TODO: here i changed the sign, check again! (victor)
    x = -Q.transpose()*result;
    
    return;
}
}}

