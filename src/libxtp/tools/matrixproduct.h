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

#ifndef _VOTCA_XTP_MATRIXPRODUCT_H
#define _VOTCA_XTP_MATRIXPRODUCT_H

#include <stdio.h>
#include <boost/format.hpp>
#include <votca/xtp/elements.h>
#include <votca/ctp/logger.h>
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <votca/tools/constants.h>

namespace votca {
    namespace xtp {

        using namespace std;
        namespace ub = boost::numeric::ublas;

        class MatProd : public ctp::QMTool {
        public:

            MatProd() {
            };

            ~MatProd() {
            };

            string Identify() {
                return "matprod";
            }

            void Initialize(Property *options);
            bool Evaluate();


        private:


            ctp::Logger _log;

        };

        void MatProd::Initialize(Property* options) {

            //UpdateWithDefaults( options, "xtp" );

            return;
        }



bool MatProd::Evaluate() {

    _log.setReportLevel( ctp::logDEBUG );

    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ...");



    ub::matrix<double> _A = ub::zero_matrix<double>(2600,2600);
    ub::matrix<double> _B = ub::zero_matrix<double>(2600,2600);

    for (int i=0; i< 2600; i++)
    {
        for (int j=0; j<2600;j++){

            _A(i,j) = sqrt(double(i)*double(j));
            _B(i,j) = cos(double(i+j));
        }
    }

    std::cout << ctp::TimeStamp() << "Start" << std::endl;
    ub::matrix<double> _product = ub::prod(_A,_B);
    std::cout << ctp::TimeStamp() << "END" << std::endl;
    return true;
        }



}}


#endif
