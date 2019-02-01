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

#ifndef VOTCA_XTP_EIGEN_H
#define	VOTCA_XTP_EIGEN_H
#include <votca/tools/eigen.h>
#include <votca/xtp/votca_config.h>
#if (GWBSE_DOUBLE)
#define real_gwbse double
#else
#define real_gwbse float
#endif


namespace votca {
    namespace xtp {
      
 typedef Eigen::Matrix<real_gwbse, Eigen::Dynamic, Eigen::Dynamic> MatrixXfd;
 typedef Eigen::Matrix<real_gwbse, Eigen::Dynamic, 1> VectorXfd;
 typedef Eigen::Matrix<double,9,1> Vector9d;

    }}





#endif	// VOTCA_XTP_EIGEN_H 

