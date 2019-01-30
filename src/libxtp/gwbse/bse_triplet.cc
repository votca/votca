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

#include <votca/xtp/orbitals.h>
#include <votca/xtp/bse_triplet.h>
#include <votca/xtp/bse_operator.h>

namespace votca {
	namespace xtp {

		BSE_Triplet::BSE_Triplet(Orbitals& orbitals,ctp::Logger &log,TCMatrix_gwbse& Mmn,const Eigen::MatrixXd& Hqp) :
			BSE_OPERATOR(orbitals, log, Mmn, Hqp) {}

		Eigen::VectorXd BSE_Triplet::col(int index) const
		{
			std::cout << std::endl << " Get Column Triplet " << index << std::endl;
			return BSE_Triplet::Hqp_col(index) + BSE_Triplet::Hd_col(index);
		}

	}
}

