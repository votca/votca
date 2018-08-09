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

#include "votca/xtp/orbitals.h"
#include <votca/xtp/version.h>
#include <votca/tools/elements.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>



using namespace votca::tools;

namespace votca {
    namespace xtp {

        Orbitals::Orbitals() {

            _basis_set_size = 0;
            _occupied_levels = 0;
            _unoccupied_levels = 0;
            _number_of_electrons = 0;
            _self_energy = 0.0;
            _qm_energy = 0.0;
            _ECP = "";
            _bsetype = "";

            // GW-BSE
            _qpmin = 0;
            _qpmax = 0;
            _qptotal = 0;

            _rpamin = 0;
            _rpamax = 0;

            _ScaHFX = 0.0;

            _bse_cmin = 0;
            _bse_cmax = 0;
            _bse_vmin = 0;
            _bse_vmax = 0;
            _bse_vtotal = 0;
            _bse_ctotal = 0;
            _bse_size = 0;
            _bse_nmax = 0;

        };

        Orbitals::~Orbitals() {
            std::vector< QMAtom* >::iterator it;
            for (it = _atoms.begin(); it != _atoms.end(); ++it) delete *it;
        };

        void Orbitals::setNumberOfLevels(int occupied_levels,int unoccupied_levels) {
            _occupied_levels = occupied_levels;
            _unoccupied_levels = unoccupied_levels;
        }
      
         /**
         *
         * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
         * @return vector with indices off all aorbitals degenerate to this including itself
         */
        std::vector<int> Orbitals::CheckDegeneracy(int level, double energy_difference)const{
          
          std::vector<int> result=std::vector<int>(0);
          if(level>_mo_energies.size()){
            throw std::runtime_error("Level for degeneracy is higher than maximum level");
          }
          double MOEnergyLevel =_mo_energies(level);

                for (int i =0; i < _mo_energies.size(); ++i) {
                    if (std::abs(_mo_energies(i) - MOEnergyLevel) * tools::conv::hrt2ev < energy_difference) {
                      result.push_back(i);
                    }
                }
          return result;
        }

        std::vector<int> Orbitals::SortEnergies() {
            if (tools::globals::verbose) cout << "... ... Sorting energies" << endl;
            std::vector<int>index = std::vector<int>(_mo_energies.size());
            std::iota(index.begin(), index.end(), 0);
            std::stable_sort(index.begin(), index.end(), [this](int i1, int i2) {
                return this->MOEnergies()[i1] > this->MOEnergies()[i2];
            });
            return index;
        }

        /// Writes a PDB file
        //TODO move to filewriter PDB

        void Orbitals::WriteXYZ(const std::string& filename, string header) const{
            
          std::ofstream out(filename);
          if (!out.is_open()) {
                throw std::runtime_error("Bad file handle: " + filename);
            }
          out<<_atoms.size()<<endl;
          out<<header<<endl;
          for (const QMAtom* atom:_atoms) {
                const tools::vec pos = atom->getPos() * tools::conv::bohr2ang;
                out<<atom->getType()<<" "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<endl;
          }
          out.close();
          return;
        }

        // Determine ground state density matrix

        Eigen::MatrixXd Orbitals::DensityMatrixGroundState() const{

            Eigen::MatrixXd occstates = _mo_coefficients.block(0, 0, _mo_coefficients.rows(), _occupied_levels);
            Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
            return dmatGS;
        }

        Eigen::MatrixXd Orbitals::LambdaMatrixQuasiParticle() const{
            return _QPdiag_coefficients * _mo_coefficients.block(0, _qpmin, _mo_coefficients.rows(), _qptotal);
        }

        // Determine QuasiParticle Density Matrix

        Eigen::MatrixXd Orbitals::DensityMatrixQuasiParticle(int state) const{
            Eigen::MatrixXd lambda = LambdaMatrixQuasiParticle();
            Eigen::MatrixXd dmatQP = lambda.col(state) * lambda.col(state).transpose();
            return dmatQP;
        }

        Orbitals::Index2MO Orbitals::BSEIndex2MOIndex()const{
          Index2MO result;
          for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
            for (unsigned _c = 0; _c < _bse_ctotal; _c++) {
              result.I2v.push_back(_bse_vmin + _v);
              result.I2c.push_back(_bse_cmin + _c);
            }
          }
          return result;
        }

        Eigen::MatrixXd Orbitals::TransitionDensityMatrix(const string& spin, int state) const{
            if (!(spin == "singlet" || spin == "triplet")) {
                throw runtime_error("Spin type not known for density matrix. Available are singlet and triplet");
            }
            const MatrixXfd& _BSECoefs = (spin == "singlet") ? _BSE_singlet_coefficients : _BSE_triplet_coefficients;
            if(_BSECoefs.cols()<state || _BSECoefs.rows()<2){
                throw runtime_error("Orbitals object has no information about that state");
            }
            Eigen::MatrixXd dmatTS = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
            double sqrt2 = sqrt(2.0);
            /*Trying to implement D_{alpha,beta}= sqrt2*sum_{i}^{occ}sum_{j}^{virt}{BSEcoef(i,j)*MOcoef(alpha,i)*MOcoef(beta,j)} */
            // c stands for conduction band and thus virtual orbitals
            // v stand for valence band and thus occupied orbitals
            
            // indexing info BSE vector index to occupied/virtual orbital
            Index2MO index=BSEIndex2MOIndex();

            if (_bsetype == "full" && spin == "singlet") {
                const MatrixXfd& _BSECoefs_AR = _BSE_singlet_coefficients_AR;
#pragma omp parallel for
                for (unsigned a = 0; a < dmatTS.rows(); a++) {
                    for (unsigned b = 0; b < dmatTS.cols(); b++) {
                        for (unsigned i = 0; i < _bse_size; i++) {
                            int occ = index.I2v[i];
                            int virt =index.I2c[i];
                            dmatTS(a, b) += sqrt2 * (_BSECoefs(i, state) + _BSECoefs_AR(i, state)) * _mo_coefficients(a, occ) * _mo_coefficients(b, virt); //check factor 2??
                        }
                    }
                }
            } else {

#pragma omp parallel for
                for (unsigned a = 0; a < dmatTS.rows(); a++) {
                    for (unsigned b = 0; b < dmatTS.cols(); b++) {
                        for (unsigned i = 0; i < _bse_size; i++) {
                            int occ = index.I2v[i];
                            int virt =index.I2c[i];
                            dmatTS(a, b) += sqrt2 * _BSECoefs(i, state) * _mo_coefficients(a, occ) * _mo_coefficients(b, virt); //check factor 2??
                        }
                    }
                }

            }
            return dmatTS;
        }

        std::vector<Eigen::MatrixXd > Orbitals::DensityMatrixExcitedState(const string& spin, int state) const{
            std::vector<Eigen::MatrixXd > dmat = DensityMatrixExcitedState_R(spin, state);
            if (_bsetype == "full" && spin == "singlet") {
                std::vector<Eigen::MatrixXd > dmat_AR = DensityMatrixExcitedState_AR(spin, state);
                dmat[0] -= dmat_AR[0];
                dmat[1] -= dmat_AR[1];
            }
            return dmat;
        }

        // Excited state density matrix

        std::vector<Eigen::MatrixXd> Orbitals::DensityMatrixExcitedState_R(const string& spin, int state) const{
            if (!(spin == "singlet" || spin == "triplet")) {
                throw runtime_error("Spin type not known for density matrix. Available are singlet and triplet");
            }

            const MatrixXfd & BSECoefs = (spin == "singlet") ? _BSE_singlet_coefficients : _BSE_triplet_coefficients;
            if(BSECoefs.cols()<state || BSECoefs.rows()<2){
                throw runtime_error("Orbitals object has no information about that state");
            }
            /******
             *
             *    Density matrix for GW-BSE based excitations
             *
             *    - electron contribution
             *      D_ab = \sum{vc} \sum{c'} A_{vc}A_{vc'} mo_a(c)mo_b(c')
             *
             *    - hole contribution
             *      D_ab = \sum{vc} \sum{v'} A_{vc}A_{v'c} mo_a(v)mo_b(v')
             *
             *
             *   more efficient:
             *
             *   - electron contribution
             *      D_ab = \sum{c} \sum{c'} mo_a(c)mo_b(c') [ \sum{v} A_{vc}A_{vc'} ]
             *           = \sum{c} \sum{c'} mo_a(c)mo_b(c') A_{cc'}
             *
             *   - hole contribution
             *      D_ab = \sum{v} \sum{v'} mo_a(v)mo_b(v') [ \sum{c} A_{vc}A_{v'c} ]
             *           = \sum{v} \sum{v'} mo_a(v)mo_b(v') A_{vv'}
             *
             */
            std::vector<Eigen::MatrixXd > dmatEX;
            dmatEX.resize(2);
            dmatEX[0] = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            dmatEX[1] = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            Index2MO index=BSEIndex2MOIndex();

            // electron assist matrix A_{cc'}
            Eigen::MatrixXd Acc = Eigen::MatrixXd::Zero(_bse_ctotal, _bse_ctotal);
            Eigen::MatrixXd Avv = Eigen::MatrixXd::Zero(_bse_vtotal, _bse_vtotal);

            for (unsigned idx1 = 0; idx1 < _bse_size; idx1++) {
                int v = index.I2v[idx1];
                int c = index.I2c[idx1];
                // electron assist matrix A_{cc'}
#pragma omp parallel for
                for (unsigned _c2 = _bse_cmin; _c2 <= _bse_cmax; _c2++) {
                    unsigned _idx2 = (_bse_cmax - _bse_cmin + 1)*(v - _bse_vmin)+(_c2 - _bse_cmin);
                    Acc(c - _bse_cmin, _c2 - _bse_cmin) += BSECoefs(idx1, state) * BSECoefs(_idx2, state);
                }

                // hole assist matrix A_{vv'}
#pragma omp parallel for
                for (unsigned v2 = _bse_vmin; v2 <= _bse_vmax; v2++) {
                    unsigned idx2 = (_bse_cmax - _bse_cmin + 1)*(v2 - _bse_vmin)+(c - _bse_cmin);
                    Avv(v - _bse_vmin, v2 - _bse_vmin) += BSECoefs(idx1, state) * BSECoefs(idx2, state);
                }
            }

            // hole part as matrix products
            Eigen::MatrixXd occlevels = _mo_coefficients.block(0, _bse_vmin, _mo_coefficients.rows(), _bse_vtotal);
            dmatEX[0] = occlevels * Avv * occlevels.transpose();

            // electron part as matrix products
            Eigen::MatrixXd virtlevels = _mo_coefficients.block(0, _bse_cmin, _mo_coefficients.rows(), _bse_ctotal);
            dmatEX[1] = virtlevels * Acc * virtlevels.transpose();
            return dmatEX;
        }

        // Excited state density matrix

        std::vector<Eigen::MatrixXd > Orbitals::DensityMatrixExcitedState_AR(const string& spin, int state) const{
            if (!(spin == "singlet")) {
                throw runtime_error("Spin type not known for density matrix. Available is singlet");
            }
            
            const MatrixXfd& BSECoefs_AR = _BSE_singlet_coefficients_AR;
            if(BSECoefs_AR.cols()<state || BSECoefs_AR.rows()<2){
                throw runtime_error("Orbitals object has no information about that state");
            }
            /******
             *
             *    Density matrix for GW-BSE based excitations
             *
             *    - electron contribution
             *      D_ab = \sum{vc} \sum{v'} B_{vc}B_{v'c} mo_a(v)mo_b(v')
             *
             *    - hole contribution
             *      D_ab = \sum{vc} \sum{c'} B_{vc}B_{vc'} mo_a(c)mo_b(c')
             *
             *
             *   more efficient:
             *
             *   - electron contribution
             *      D_ab = \sum{v} \sum{v'} mo_a(v)mo_b(v') [ \sum{c} B_{vc}B_{v'c} ]
             *           = \sum{v} \sum{v'} mo_a(v)mo_b(v') B_{vv'}
             *
             *   - hole contribution
             *      D_ab = \sum{c} \sum{c'} mo_a(c)mo_b(c') [ \sum{v} B_{vc}B_{vc'} ]
             *           = \sum{c} \sum{c'} mo_a(c)mo_b(c') B_{cc'}
             *
             */
            std::vector<Eigen::MatrixXd > dmatAR;
            dmatAR.resize(2);
            dmatAR[0] = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            dmatAR[1] = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            Index2MO index=BSEIndex2MOIndex();

            // hole assist matrix B_{cc'}
            Eigen::MatrixXd Bcc = Eigen::MatrixXd::Zero(_bse_ctotal, _bse_ctotal);
            Eigen::MatrixXd Bvv = Eigen::MatrixXd::Zero(_bse_vtotal, _bse_vtotal);

            for (unsigned idx1 = 0; idx1 < _bse_size; idx1++) {
                int v = index.I2v[idx1];
                int c = index.I2c[idx1];
                // hole assist matrix B_{cc'}
#pragma omp parallel for
                for (unsigned c2 = _bse_cmin; c2 <= _bse_cmax; c2++) {
                    unsigned idx2 = (_bse_cmax - _bse_cmin + 1)*(v - _bse_vmin)+(c2 - _bse_cmin);
                    Bcc(c - _bse_cmin, c2 - _bse_cmin) += BSECoefs_AR(idx1, state) * BSECoefs_AR(idx2, state);
                }

                // electron assist matrix B_{vv'}
#pragma omp parallel for
                for (unsigned v2 = _bse_vmin; v2 <= _bse_vmax; v2++) {
                    unsigned idx2 = (_bse_cmax - _bse_cmin + 1)*(v2 - _bse_vmin)+(c - _bse_cmin);
                    Bvv(v - _bse_vmin, v2 - _bse_vmin) += BSECoefs_AR(idx1, state) * BSECoefs_AR(idx2, state);
                }
            }

            // hole part as matrix products
            Eigen::MatrixXd occlevels = _mo_coefficients.block(0, _bse_vmin, _mo_coefficients.rows(), _bse_vtotal);
            dmatAR[0] = occlevels * Bvv * occlevels.transpose();
            // electron part as matrix products
            Eigen::MatrixXd virtlevels = _mo_coefficients.block(0, _bse_cmin, _mo_coefficients.rows(), _bse_ctotal);
            dmatAR[1] = virtlevels * Bcc * virtlevels.transpose();
            return dmatAR;
        }

        Eigen::VectorXd Orbitals::LoewdinPopulation(const Eigen::MatrixXd & densitymatrix, const Eigen::MatrixXd & overlapmatrix, int frag){

            Eigen::VectorXd fragmentCharges = Eigen::VectorXd::Zero(2);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
            es.compute(overlapmatrix);
            Eigen::MatrixXd sqrtm1 = es.operatorInverseSqrt();
            Eigen::MatrixXd prodmat = sqrtm1 * densitymatrix*sqrtm1;

            for (int _i = 0; _i < frag; _i++) {
                fragmentCharges(0) += prodmat(_i, _i);
            }
            for (unsigned _i = frag; _i < overlapmatrix.rows(); _i++) {
                fragmentCharges(1) += prodmat(_i, _i);
            }

            return fragmentCharges;
        }

        std::vector<double> Orbitals::Oscillatorstrengths() const{
            std::vector<double> oscs;
            unsigned size = _transition_dipoles.size();
            if (size > _BSE_singlet_energies.size()) {
                size = _BSE_singlet_energies.size();
            }
            for (unsigned i = 0; i < size; ++i) {
                double osc = (_transition_dipoles[i] * _transition_dipoles[i]) * 2.0 / 3.0 * (_BSE_singlet_energies(i));
                oscs.push_back(osc);
            }
            return oscs;
        }

        double Orbitals::getTotalExcitedStateEnergy(const string& spintype, int opt_state) const{

            // total energy of the excited state
            double total_energy;
            double omega = 0.0;

            double dft_energy = getQMEnergy();

            if (spintype == "singlet") {
              if(BSESingletEnergies().size()<opt_state){
                throw std::runtime_error("Orbitals::getTotalEnergy You want a singlet which has not been calculated");
              }
                omega = BSESingletEnergies()[opt_state - 1];
            } else if (spintype == "triplet") {
               if(BSETripletEnergies().size()<opt_state){
                throw std::runtime_error("Orbitals::getTotalEnergy You want a triplet which has not been calculated");
              }
                omega = BSETripletEnergies()[opt_state - 1];
            } else {
                throw std::runtime_error("GetTotalEnergy only knows spintypes:singlet,triplet");
            }

            // DFT total energy is stored in eV
            // singlet energies are stored in Hrt...
            return total_energy = dft_energy * tools::conv::ev2hrt + omega; //  e.g. hartree
        }

        Eigen::VectorXd Orbitals::FragmentNuclearCharges(int frag) const{
         
            if (frag < 0) {
                throw runtime_error("Orbitals::FragmentNuclearCharges Fragment index is smaller than zero");
            }

            Eigen::VectorXd fragmentNuclearCharges = Eigen::VectorXd::Zero(2);
            int id = 0;
            for (const QMAtom* atom :_atoms) {
                id++;
                // get element type and determine its nuclear charge
                double crg = atom->getNuccharge();
                // add to either fragment
                if (id <= frag) {
                    fragmentNuclearCharges(0) += crg;
                } else {
                    fragmentNuclearCharges(1) += crg;
                }
            }
            return fragmentNuclearCharges;
        }

                /**
         * \brief Guess for a dimer based on monomer orbitals
         *
         * Given two monomer orbitals (A and B) constructs a guess for dimer
         * orbitals: | A 0 | and energies: [EA, EB]
         *           | 0 B |
         */
        void Orbitals::PrepareDimerGuess(const Orbitals& orbitalsA,const Orbitals& orbitalsB, Orbitals& orbitalsAB) {

            // constructing the direct product orbA x orbB
            int basisA = orbitalsA.getBasisSetSize();
            int basisB = orbitalsB.getBasisSetSize();

            int levelsA = orbitalsA.getNumberOfLevels();
            int levelsB = orbitalsB.getNumberOfLevels();

            int electronsA = orbitalsA.getNumberOfElectrons();
            int electronsB = orbitalsB.getNumberOfElectrons();

            Eigen::MatrixXd& mo_coefficients = orbitalsAB.MOCoefficients();
            mo_coefficients = Eigen::MatrixXd::Zero(basisA + basisB, levelsA + levelsB);

            // AxB = | A 0 |  //   A = [EA, EB]  //
            //       | 0 B |  //                 //
            orbitalsAB.setBasisSetSize(basisA + basisB);
            orbitalsAB.setNumberOfLevels(electronsA - electronsB,
                    levelsA + levelsB - electronsA - electronsB);
            orbitalsAB.setNumberOfElectrons(electronsA + electronsB);

            mo_coefficients.block(0, 0, basisA, levelsA) = orbitalsA.MOCoefficients();
            mo_coefficients.block(basisA, levelsA, basisB, levelsB) = orbitalsB.MOCoefficients();

            Eigen::VectorXd& energies = orbitalsAB.MOEnergies();
            energies.resize(levelsA + levelsB);

            energies.segment(0, levelsA) = orbitalsA.MOEnergies();
            energies.segment(levelsA, levelsB) = orbitalsB.MOEnergies();
            return;
        }
        //TODO move to Filereader

        void Orbitals::LoadFromXYZ(const std::string& filename) {

            string line;
            std::ifstream in;
            string type;
            in.open(filename.c_str(), std::ios::in);
            if (!in) throw runtime_error(string("Error reading coordinates from: ")
                    + filename);
            int atomCount = 0;
            std::getline(in, line);
            
            Tokenizer tok1(line," \t");
            std::vector<std::string> line1;
            tok1.ToVector(line1);
            if(line1.size()!=1){
              throw std::runtime_error("First line of xyz file should contain number of atoms, nothing else.");
            }
            std::getline(in, line);//Comment line
            
            if (in.is_open()) {
                while (in.good()) {
                    std::getline(in, line);

                    vector< string > split;
                    Tokenizer toker(line, " \t");
                    toker.ToVector(split);
                    if (!split.size() ||
                            split.size() != 4 ||
                            split[0] == "#" ||
                            split[0].substr(0, 1) == "#") {
                        continue;
                    }
                    // Interesting information written here: e.g. 'C 0.000 0.000 0.000'

                    string element = split[0];
                    double x = boost::lexical_cast<double>(split[1]);
                    double y = boost::lexical_cast<double>(split[2]);
                    double z = boost::lexical_cast<double>(split[3]);
                    tools::vec pos = tools::vec(x, y, z);
                    AddAtom(atomCount, element, pos * tools::conv::ang2bohr);
                    atomCount++;

                }
            } else {
                throw std::runtime_error("No such file: '" + filename + "'.");
            }
            return;
        }

        void Orbitals::WriteToCpt(const std::string& filename) const{
            CheckpointFile cpf(filename, true);
            WriteToCpt(cpf);
        }

        void Orbitals::WriteToCpt(CheckpointFile f) const{
            std::string name = "QMdata";
            CptLoc orbGr = f.getHandle().createGroup("/" + name);
            WriteToCpt(orbGr);
        }

        void Orbitals::WriteToCpt(CptLoc parent) const{
            try {
                CheckpointWriter w(parent);
                w(XtpVersionStr(), "Version");
                w(_basis_set_size, "basis_set_size");
                w(_occupied_levels, "occupied_levels");
                w(_unoccupied_levels, "unoccupied_levels");
                w(_number_of_electrons, "number_of_electrons");

                w(_mo_energies, "mo_energies");
                w(_mo_coefficients, "mo_coefficients");

                // write qmatoms

                {
                    CptLoc qmAtomsGr = parent.createGroup("qmatoms");
                    size_t count = 0;
                    for (const auto& qma : _atoms) {
                        CptLoc tempLoc = qmAtomsGr.createGroup("atom" + std::to_string(count));
                        qma->WriteToCpt(tempLoc);
                        ++count;
                    }

                }

                w(_qm_energy, "qm_energy");
                w(_qm_package, "qm_package");
                w(_self_energy, "self_energy");

                w(_dftbasis, "dftbasis");
                w(_auxbasis, "auxbasis");

                w(_rpamin, "rpamin");
                w(_rpamax, "rpamax");
                w(_qpmin, "qpmin");
                w(_qpmax, "qpmax");
                w(_bse_vmin, "bse_vmin");
                w(_bse_vmax, "bse_vmax");
                w(_bse_cmin, "bse_cmin");
                w(_bse_cmax, "bse_cmax");

                w(_ScaHFX, "ScaHFX");

                w(_bsetype, "bsetype");
                w(_ECP, "ECP");

                w(_QPpert_energies, "QPpert_energies");
                w(_QPdiag_energies, "QPdiag_energies");

                w(_QPdiag_coefficients, "QPdiag_coefficients");
                w(_eh_t, "eh_t");

                w(_eh_s, "eh_s");

                w(_BSE_singlet_energies, "BSE_singlet_energies");

                w(_BSE_singlet_coefficients, "BSE_singlet_coefficients");

                w(_BSE_singlet_coefficients_AR, "BSE_singlet_coefficients_AR");

                w(_transition_dipoles, "transition_dipoles");

                w(_BSE_triplet_energies, "BSE_triplet_energies");
                w(_BSE_triplet_coefficients, "BSE_triplet_coefficients");

            } catch (H5::Exception& error) {
                throw std::runtime_error(error.getDetailMsg());
            }

        }

        void Orbitals::ReadFromCpt(const std::string& filename) {
            CheckpointFile cpf(filename, false);
            ReadFromCpt(cpf);
        }

        void Orbitals::ReadFromCpt(CheckpointFile f) {
            std::string name = "QMdata";
            CptLoc orbGr = f.getHandle().openGroup("/" + name);
            ReadFromCpt(orbGr);
        }

        void Orbitals::ReadFromCpt(CptLoc parent) {
            try {
                CheckpointReader r(parent);
                r(_basis_set_size, "basis_set_size");
                r(_occupied_levels, "occupied_levels");
                r(_unoccupied_levels, "unoccupied_levels");
                r(_number_of_electrons, "number_of_electrons");

                r(_mo_energies, "mo_energies");
                r(_mo_coefficients, "mo_coefficients");
                // Read qmatoms
                {
                    CptLoc qmAtomsGr = parent.openGroup("qmatoms");
                    size_t count = qmAtomsGr.getNumObjs();
                    if(this->QMAtoms().size()>0){
                      std::vector< QMAtom* >::iterator it;
                      for (it = _atoms.begin(); it != _atoms.end(); ++it) delete *it;
                      _atoms.clear();
                    }

                    for (size_t i = 0; i < count; ++i) {
                        CptLoc tempLoc = qmAtomsGr.openGroup("atom" + std::to_string(i));
                        QMAtom temp;
                        temp.ReadFromCpt(tempLoc);
                        AddAtom(temp);
                    }
                }

                r(_qm_energy, "qm_energy");
                r(_qm_package, "qm_package");
                r(_self_energy, "self_energy");

                r(_dftbasis, "dftbasis");
                r(_auxbasis, "auxbasis");

                r(_rpamin, "rpamin");
                r(_rpamax, "rpamax");
                r(_qpmin, "qpmin");
                r(_qpmax, "qpmax");
                r(_bse_vmin, "bse_vmin");
                r(_bse_vmax, "bse_vmax");
                r(_bse_cmin, "bse_cmin");
                r(_bse_cmax, "bse_cmax");
                _bse_vtotal = _bse_vmax - _bse_vmin + 1;
                _bse_ctotal = _bse_cmax - _bse_cmin + 1;
                _bse_size = _bse_vtotal * _bse_ctotal;

                r(_ScaHFX, "ScaHFX");

                r(_bsetype, "bsetype");
                r(_ECP, "ECP");

                r(_QPpert_energies, "QPpert_energies");
                r(_QPdiag_energies, "QPdiag_energies");

                r(_QPdiag_coefficients, "QPdiag_coefficients");
                r(_eh_t, "eh_t");

                r(_eh_s, "eh_s");

                r(_BSE_singlet_energies, "BSE_singlet_energies");

                r(_BSE_singlet_coefficients, "BSE_singlet_coefficients");

                r(_BSE_singlet_coefficients_AR, "BSE_singlet_coefficients_AR");

                r(_transition_dipoles, "transition_dipoles");

                r(_BSE_triplet_energies, "BSE_triplet_energies");
                r(_BSE_triplet_coefficients, "BSE_triplet_coefficients");

            } catch (H5::Exception& error) {
                throw std::runtime_error(error.getDetailMsg());
            }
        }
    }
}
