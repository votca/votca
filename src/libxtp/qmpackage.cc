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

// Overload of uBLAS prod function with MKL/GSL implementations
#include "votca/xtp/qmpackage.h"
#include "votca/xtp/aomatrix.h"
#include <votca/xtp/molecule.h>
#include <votca/xtp/atom.h>
#include <boost/algorithm/string.hpp>

namespace votca {
    namespace xtp {
      using std::flush;
        void QMPackage::ReorderOutput(Orbitals& orbitals) {
            BasisSet dftbasisset;
            dftbasisset.LoadBasisSet(_basisset_name);
            if (!orbitals.hasQMAtoms()) {
                throw std::runtime_error("Orbitals object has no QMAtoms");
            }

            AOBasis dftbasis;
            dftbasis.AOBasisFill(dftbasisset, orbitals.QMAtoms());
            //necessary to update nuclear charges on qmatoms
            if (_write_pseudopotentials) {
                BasisSet ecps;
                ecps.LoadPseudopotentialSet(_ecp_name);
                AOBasis ecpbasis;
                ecpbasis.ECPFill(ecps, orbitals.QMAtoms());
            }

            if (orbitals.hasAOOverlap()) {
                dftbasis.ReorderMatrix(orbitals.AOOverlap(), getPackageName(), "xtp");
                XTP_LOG(logDEBUG, *_pLog) << "Reordered Overlap matrix" << flush;
            }
            if (orbitals.hasAOVxc()) {
                dftbasis.ReorderMatrix(orbitals.AOVxc(), getPackageName(), "xtp");
                XTP_LOG(logDEBUG, *_pLog) << "Reordered VXC matrix" << flush;
            }
            if (orbitals.hasMOCoefficients()) {
                dftbasis.ReorderMOs(orbitals.MOCoefficients(), getPackageName(), "xtp");
                XTP_LOG(logDEBUG, *_pLog) << "Reordered MOs" << flush;
            }

            return;
        }

        void QMPackage::ReorderMOsBack(Orbitals& orbitals) {
            BasisSet dftbasisset;
            dftbasisset.LoadBasisSet(_basisset_name);
            if (!orbitals.hasQMAtoms()) {
                throw std::runtime_error("Orbitals object has no QMAtoms");
            }
            AOBasis dftbasis;
            dftbasis.AOBasisFill(dftbasisset, orbitals.QMAtoms());
            dftbasis.ReorderMOs(orbitals.MOCoefficients(), "xtp", getPackageName());
            return;
        }

        std::vector<QMPackage::MinimalMMCharge > QMPackage::SplitMultipoles(const PolarSite& aps) {

            std::vector< QMPackage::MinimalMMCharge > multipoles_split;
            // Calculate virtual charge positions
            double a = _dpl_spacing; // this is in a0
            double mag_d = aps.getDipole().norm();// this is in e * a0
            const Eigen::Vector3d dir_d = aps.getDipole().normalized();
            const Eigen::Vector3d A = aps.getPos() + 0.5 * a * dir_d; // converted to AA
            const Eigen::Vector3d B = aps.getPos() - 0.5 * a * dir_d;
            double qA = mag_d / a;
            double qB = -qA;
            multipoles_split.push_back(MinimalMMCharge(A, qA));
            multipoles_split.push_back(MinimalMMCharge(B, qB));


            if (aps.getRank() > 1) {
                const Eigen::Matrix3d components = aps.CalculateCartesianMultipole();
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
                es.computeDirect(components);
                double a = 2 * _dpl_spacing;
                for (int i = 0; i < 3; i++) {
                    double q = es.eigenvalues()[i] / (a * a);
                    const Eigen::Vector3d vec1 = aps.getPos() + 0.5 * a * es.eigenvectors().col(i);
                    const Eigen::Vector3d vec2 = aps.getPos() - 0.5 * a * es.eigenvectors().col(i);
                    multipoles_split.push_back(MinimalMMCharge(vec1, q));
                    multipoles_split.push_back(MinimalMMCharge(vec2, q));
                }
            }
            return multipoles_split;
        }
        
      void QMPackage::setMultipoleBackground(const std::shared_ptr<MMRegion>& PolarSegments ) {
        if(PolarSegments->size()==0){
          std::cout<<"WARNING::The Multipole Background has no entries!"<<std::endl;
          return;
        }
      _PolarSegments = PolarSegments;
      _write_charges = true;
      
      WriteChargeOption();
    }
      
      std::vector<std::string> QMPackage::GetLineAndSplit(std::ifstream& input_file,const std::string separators ){
          std::string line;
          getline(input_file, line);
          boost::trim(line);
          std::vector<std::string> row;
          boost::algorithm::split(row, line, boost::is_any_of(separators), boost::algorithm::token_compress_on);
          return row;
        }


    }
}
