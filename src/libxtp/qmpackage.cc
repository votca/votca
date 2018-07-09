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

namespace votca {
    namespace xtp {
      using std::flush;
        void QMPackage::ReorderOutput(Orbitals* _orbitals) {
            BasisSet _dftbasisset;
            _dftbasisset.LoadBasisSet(_basisset_name);
            if (!_orbitals->hasQMAtoms()) {
                throw std::runtime_error("Orbitals object has no QMAtoms");
            }


            AOBasis _dftbasis;
            _dftbasis.AOBasisFill(&_dftbasisset, _orbitals->QMAtoms());
            //necessary to update nuclear charges on qmatoms
            if (_write_pseudopotentials) {
                BasisSet _ecps;
                _ecps.LoadPseudopotentialSet(_ecp_name);
                AOBasis _ecpbasis;
                _ecpbasis.ECPFill(&_ecps, _orbitals->QMAtoms());
            }

            if (_orbitals->hasAOOverlap()) {
                _dftbasis.ReorderMatrix(_orbitals->AOOverlap(), getPackageName(), "xtp");
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered Overlap matrix" << flush;
            }
            if (_orbitals->hasAOVxc()) {
                _dftbasis.ReorderMatrix(_orbitals->AOVxc(), getPackageName(), "xtp");
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered VXC matrix" << flush;
            }
            if (_orbitals->hasMOCoefficients()) {
                _dftbasis.ReorderMOs(_orbitals->MOCoefficients(), getPackageName(), "xtp");
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered MOs" << flush;
            }

            return;
        }

        void QMPackage::ReorderMOsBack(Orbitals* _orbitals) {
            BasisSet _dftbasisset;
            _dftbasisset.LoadBasisSet(_basisset_name);
            if (!_orbitals->hasQMAtoms()) {
                throw std::runtime_error("Orbitals object has no QMAtoms");
            }
            AOBasis _dftbasis;
            _dftbasis.AOBasisFill(&_dftbasisset, _orbitals->QMAtoms());
            _dftbasis.ReorderMOs(_orbitals->MOCoefficients(), "xtp", getPackageName());
            return;
        }

        std::vector<std::vector<double> > QMPackage::SplitMultipoles(ctp::APolarSite* aps) {

            std::vector< std::vector<double> > multipoles_split;

            const tools::vec pos = aps->getPos() * tools::conv::nm2ang;
            tools::vec tot_dpl = tools::vec(0.0);
            if (_with_polarization) {
                tot_dpl += aps->getU1();
            }
            if (aps->getRank() > 0) {
                tot_dpl += aps->getQ1();
            }
            // Calculate virtual charge positions
            double a = _dpl_spacing; // this is in nm
            double mag_d = abs(tot_dpl); // this is in e * nm
            if (mag_d > 1e-9) {
                tools::vec dir_d = tot_dpl.normalize();
                tools::vec A = pos + 0.5 * a * dir_d * tools::conv::nm2ang; // converted to AA
                tools::vec B = pos - 0.5 * a * dir_d * tools::conv::nm2ang;
                double qA = mag_d / a;
                double qB = -qA;
                multipoles_split.push_back({A.getX(), A.getY(), A.getZ(), qA});
                multipoles_split.push_back({B.getX(), B.getY(), B.getZ(), qB});
            }


            if (aps->getRank() > 1) {
                tools::matrix components = aps->getQ2cartesian();
                tools::matrix::eigensystem_t system;
                components.SolveEigensystem(system);
                double a = 2 * _dpl_spacing;
                for (unsigned i = 0; i < 3; i++) {

                    double q = system.eigenvalues[i] / (a * a);
                    if (std::abs(q) < 1e-9) {
                        continue;
                    }
                    tools::vec vec1 = pos + 0.5 * a * system.eigenvecs[i] * tools::conv::nm2ang;
                    tools::vec vec2 = pos - 0.5 * a * system.eigenvecs[i] * tools::conv::nm2ang;

                    multipoles_split.push_back({vec1.getX(), vec1.getY(), vec1.getZ(), q});
                    multipoles_split.push_back({vec2.getX(), vec2.getY(), vec2.getZ(), q});

                }

            }

            return multipoles_split;
        }

        void QMPackage::addLinkers(std::vector< ctp::Segment* > &segments, ctp::QMPair* pair, std::vector<std::string> linker_names) {
            ctp::Segment* seg1 = pair->Seg1();
            ctp::Segment* seg2 = pair->Seg2();
            ctp::Topology* _top = seg1->getTopology();
            std::vector<ctp::Segment*> segmentsInMolecule = _top->Segments();
            ctp::Molecule* moleculeSeg1 = seg1->getMolecule();
            ctp::Molecule* moleculeSeg2 = seg2->getMolecule();
            int moleculeIdSeg1 = moleculeSeg1-> getId();
            int moleculeIdSeg2 = moleculeSeg2-> getId();

            if (moleculeIdSeg1 == moleculeIdSeg2) {

                int idSeg1 = seg1->getId();
                int idSeg2 = seg2->getId();

                std::cout << "\n\nsegment size before addLinker: " << segments.size() << "\n";

                std::vector<ctp::Segment*>::iterator it;
                for (it = segmentsInMolecule.begin(); it != segmentsInMolecule.end(); ++it) {
                    ctp::Molecule* moleculeSegIt = (*it)->getMolecule();
                    int moleculeIdOfSegIt = moleculeSegIt-> getId();
                    int idIterator = (*it)->getId();
                    if (moleculeIdOfSegIt == moleculeIdSeg1 && idIterator != idSeg1 && idIterator != idSeg2
                            && isLinker((*it)->getName(), linker_names)) {
                        segments.push_back((*it));
                    }
                }

                std::cout << "\n\nsegment size after addLinker: " << segments.size() << "\n";

            }
            return;
        }

        bool QMPackage::isLinker(std::string name, std::vector< std::string> linker_names) {
            return (std::find(linker_names.begin(), linker_names.end(), name) != linker_names.end());
        }

        bool QMPackage::WriteInputFilePBC(ctp::QMPair* pair, Orbitals* orbitals, std::vector<std::string> linker_names) {

            //std::cout << "IDFT writes input with PBC" << std::endl;

            ctp::Segment* seg1 = pair->Seg1();
            ctp::Segment* seg2 = pair->Seg2();
            ctp::Segment* ghost = NULL;

            ctp::Topology* _top = seg1->getTopology();

            ctp::vec r1 = seg1->getPos();
            ctp::vec r2 = seg2->getPos();

            ctp::vec _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

            // Check whether pair formed across periodic boundary
            if (abs(r2 - r1 - _R) > 1e-8) {
                ghost = new ctp::Segment(seg2);
                //ghost->TranslateBy(r1 - r2 + _R); // DO NOT USE THIS METHOD !
                std::vector<ctp::Atom*>::iterator ait;
                for (ait = ghost->Atoms().begin(); ait != ghost->Atoms().end(); ++ait) {
                    (*ait)->setQMPos((*ait)->getQMPos() + r1 - r2 + _R);
                }
            }

            std::vector< ctp::Segment* > segments;
            segments.push_back(seg1);

            if (ghost) {
                segments.push_back(ghost);
            } else {
                segments.push_back(seg2);
            }


            std::cout << "Number of linker names " << linker_names.size() << std::endl;
            for (size_t i = 0; i < linker_names.size(); i++) {
                std::cout << linker_names[i] << std::endl;
            }

            addLinkers(segments, pair, linker_names);

            std::cout << "\n\nBefore writing: " << segments.size() << "\n";

            WriteInputFile(segments, orbitals);

            delete ghost;
            return true;
        }



    }
}
