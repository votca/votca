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
#include "votca/xtp/qmstate.h"
#include "votca/xtp/aomatrix.h"
#include <votca/xtp/version.h>
#include <votca/tools/elements.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/vc2index.h>
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
            _useTDA = false;

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
            for (QMAtom* atom:_atoms) delete atom;
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
            std::vector<int>index = std::vector<int>(_mo_energies.size());
            std::iota(index.begin(), index.end(), 0);
            std::stable_sort(index.begin(), index.end(), [this](int i1, int i2) {
                return this->MOEnergies()[i1] < this->MOEnergies()[i2];
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

        Eigen::MatrixXd Orbitals::DensityMatrixFull(const QMState& state) const{
          if(state.isTransition()){
            return this->TransitionDensityMatrix(state);
          }
          Eigen::MatrixXd result=this->DensityMatrixGroundState();
          if(state.Type().isExciton()){
             std::vector< Eigen::MatrixXd > DMAT = DensityMatrixExcitedState(state);
             result=result- DMAT[0] + DMAT[1]; // Ground state + hole_contribution + electron contribution
          }else if(state.Type()==QMStateType::DQPstate){
            Eigen::MatrixXd DMATQP=DensityMatrixQuasiParticle(state);
            if (state.Index() > getHomo()) {
              result+= DMATQP;
            } else {
              result-= DMATQP;
            }
          }else if(state.Type()!=QMStateType::Gstate){
            throw std::runtime_error("DensityMatrixFull does not yet implement QMStateType:"+state.Type().ToLongString());
          }
          return result;
        }
        
        // Determine ground state density matrix

        Eigen::MatrixXd Orbitals::DensityMatrixGroundState() const{
            if (!hasMOCoefficients()) {
                throw std::runtime_error("Orbitals file does not contain MO coefficients");
            }
            Eigen::MatrixXd occstates = _mo_coefficients.block(0, 0, _mo_coefficients.rows(), _occupied_levels);
            Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
            return dmatGS;
        }

        Eigen::MatrixXd Orbitals::CalculateQParticleAORepresentation() const{
          if (!hasQPdiag()) {
                throw std::runtime_error("Orbitals file does not contain QP coefficients");
            }
            return _mo_coefficients.block(0, _qpmin, _mo_coefficients.rows(), _qptotal)*_QPdiag_coefficients;
        }

        // Determine QuasiParticle Density Matrix

        Eigen::MatrixXd Orbitals::DensityMatrixQuasiParticle(const QMState& state) const{
          if(state.Type()!=QMStateType::PQPstate){
            throw std::runtime_error("State:"+state.ToString()+" is not a quasiparticle state");
          }
            Eigen::MatrixXd lambda = CalculateQParticleAORepresentation();
            Eigen::MatrixXd dmatQP = lambda.col(state.Index()-_qpmin) * lambda.col(state.Index()-_qpmin).transpose();
            return dmatQP;
        }
        
        Eigen::Vector3d Orbitals::CalcCoM()const{
          tools::Elements elements;
          
          Eigen::Vector3d CoM = Eigen::Vector3d::Zero();
          double totalmass=0.0;
          for (QMAtom* atom : _atoms) {
            double mass=elements.getMass(atom->getType());
            totalmass+=mass;
            CoM += mass* atom->getPos().toEigen();
          }
          CoM /= totalmass;
          return CoM;
        }
        
        Eigen::Vector3d Orbitals::CalcElDipole(const QMState& state) {
          Eigen::Vector3d CoM=CalcCoM();
          Eigen::Vector3d nuclei_dip = Eigen::Vector3d::Zero();
          if (!state.isTransition()) {
            for (QMAtom* atom : _atoms) {
              nuclei_dip += (atom->getPos().toEigen() - CoM) * atom->getNuccharge();
            }
          }

          BasisSet basis;
          basis.LoadBasisSet(this->getDFTbasis());
          AOBasis aobasis;
          aobasis.AOBasisFill(basis, _atoms);
          AODipole dipole;
          dipole.setCenter(CoM);
          dipole.Fill(aobasis);
          Eigen::MatrixXd dmat = this->DensityMatrixFull(state);
          Eigen::Vector3d electronic_dip;
          for (int i = 0; i < 3; ++i) {
            electronic_dip(i) = dmat.cwiseProduct(dipole.Matrix()[i]).sum();
          }
          return nuclei_dip - electronic_dip;
        }

        Eigen::MatrixXd Orbitals::TransitionDensityMatrix(const QMState& state) const{
            if (state.Type() != QMStateType::Singlet) {
                throw runtime_error("Spin type not known for transition density matrix. Available only for singlet");
            }
            const MatrixXfd& BSECoefs = _BSE_singlet_coefficients;
            if(BSECoefs.cols()<state.Index()+1 || BSECoefs.rows()<2){
                throw runtime_error("Orbitals object has no information about state:"+state.ToString());
            }
            
            // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
            
            /*Trying to implement D_{alpha,beta}= sqrt2*sum_{i}^{occ}sum_{j}^{virt}{BSEcoef(i,j)*MOcoef(alpha,i)*MOcoef(beta,j)} */
            // c stands for conduction band and thus virtual orbitals
            // v stand for valence band and thus occupied orbitals
#if (GWBSE_DOUBLE)
            Eigen::VectorXd coeffs= BSECoefs.col(state.Index());
#else
            Eigen::VectorXd coeffs= BSECoefs.col(state.Index()).cast<double>();
#endif
            if(!_useTDA){
#if (GWBSE_DOUBLE)
                coeffs+=_BSE_singlet_coefficients_AR.col(state.Index());
#else
            coeffs+=_BSE_singlet_coefficients_AR.col(state.Index()).cast<double>();
#endif
            }
            coeffs*=std::sqrt(2.0);
            vc2index index=vc2index(_bse_vmin,_bse_cmin,_bse_ctotal);
            Eigen::MatrixXd dmatTS = Eigen::MatrixXd::Zero(_basis_set_size, _basis_set_size);
            
            for (int i = 0; i < _bse_size; i++) {
                dmatTS.noalias()+= coeffs(i) * _mo_coefficients.col(index.v(i)) * _mo_coefficients.col(index.c(i)).transpose();
            }

            return dmatTS;
        }

        std::vector<Eigen::MatrixXd > Orbitals::DensityMatrixExcitedState(const QMState& state) const{
            std::vector<Eigen::MatrixXd > dmat = DensityMatrixExcitedState_R(state);
            if (!_useTDA && state.Type() == QMStateType::Singlet) {
                std::vector<Eigen::MatrixXd > dmat_AR = DensityMatrixExcitedState_AR(state);
                dmat[0] -= dmat_AR[0];
                dmat[1] -= dmat_AR[1];
            }
            return dmat;
        }

        // Excited state density matrix

        std::vector<Eigen::MatrixXd> Orbitals::DensityMatrixExcitedState_R(const QMState& state) const{
            if (!state.Type().isExciton()) {
                throw runtime_error("Spin type not known for density matrix. Available are singlet and triplet");
            }

            const MatrixXfd & BSECoefs = (state.Type() == QMStateType::Singlet) ? _BSE_singlet_coefficients : _BSE_triplet_coefficients;
            if(BSECoefs.cols()<state.Index()+1 || BSECoefs.rows()<2){
                throw runtime_error("Orbitals object has no information about state:"+state.ToString());
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

#if (GWBSE_DOUBLE)
            Eigen::VectorXd coeffs= BSECoefs.col(state.Index());
#else
            Eigen::VectorXd coeffs= BSECoefs.col(state.Index()).cast<double>();
#endif
         
            std::vector<Eigen::MatrixXd > dmatEX(2);
            // hole part as matrix products
            Eigen::MatrixXd occlevels = _mo_coefficients.block(0, _bse_vmin, _mo_coefficients.rows(), _bse_vtotal);
            dmatEX[0] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

            // electron part as matrix products
            Eigen::MatrixXd virtlevels = _mo_coefficients.block(0, _bse_cmin, _mo_coefficients.rows(), _bse_ctotal);
            dmatEX[1] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();

            return dmatEX;
        }

   
Eigen::MatrixXd Orbitals::CalcAuxMat_vv(const Eigen::VectorXd& coeffs)const{
    Eigen::MatrixXd Mvv = Eigen::MatrixXd::Zero(_bse_vtotal, _bse_vtotal);
    vc2index index = vc2index(_bse_vmin, _bse_cmin, _bse_ctotal);
    for (int idx1 = 0; idx1 < _bse_size; idx1++) {
        int v = index.v(idx1) - _bse_vmin;
        int c = index.c(idx1) - _bse_cmin;
#pragma omp parallel for
        for (int v2 = 0; v2 < _bse_vtotal; v2++) {
            int idx2 = index.I(v2+_bse_vmin, c+_bse_cmin);
            Mvv(v, v2) += coeffs(idx1) * coeffs(idx2);
        }
    }
    return Mvv;
}

Eigen::MatrixXd Orbitals::CalcAuxMat_cc(const Eigen::VectorXd& coeffs)const{
    Eigen::MatrixXd Mcc = Eigen::MatrixXd::Zero(_bse_ctotal, _bse_ctotal);
    vc2index index = vc2index(_bse_vmin, _bse_cmin, _bse_ctotal);
    for (int idx1 = 0; idx1 < _bse_size; idx1++) {
        int v = index.v(idx1) - _bse_vmin;
        int c = index.c(idx1) - _bse_cmin;
#pragma omp parallel for
        for (int c2 = 0; c2 < _bse_ctotal; c2++) {
            int idx2 = index.I(v+_bse_vmin, c2+_bse_cmin);
            Mcc(c, c2) += coeffs(idx1) * coeffs(idx2);
        }

    }
    return Mcc;
}

        std::vector<Eigen::MatrixXd > Orbitals::DensityMatrixExcitedState_AR(const QMState& state) const{
            if (state.Type() != QMStateType::Singlet) {
                throw runtime_error("Spin type not known for density matrix. Available is singlet");
            }
            
            const MatrixXfd& BSECoefs_AR = _BSE_singlet_coefficients_AR;
            if(BSECoefs_AR.cols()<state.Index()+1 || BSECoefs_AR.rows()<2){
                throw runtime_error("Orbitals object has no information about state:"+state.ToString());
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
            
#if (GWBSE_DOUBLE)
            Eigen::VectorXd coeffs= BSECoefs_AR.col(state.Index());
#else
            Eigen::VectorXd coeffs= BSECoefs_AR.col(state.Index()).cast<double>();
#endif

            std::vector<Eigen::MatrixXd > dmatAR(2);
            Eigen::MatrixXd virtlevels = _mo_coefficients.block(0, _bse_cmin, _mo_coefficients.rows(), _bse_ctotal);
            dmatAR[0] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();
            // electron part as matrix products
            Eigen::MatrixXd occlevels = _mo_coefficients.block(0, _bse_vmin, _mo_coefficients.rows(), _bse_vtotal);
            dmatAR[1] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

            
            return dmatAR;
        }

        Eigen::VectorXd Orbitals::LoewdinPopulation(const Eigen::MatrixXd & densitymatrix, const Eigen::MatrixXd & overlapmatrix, int frag){

            Eigen::VectorXd fragmentCharges = Eigen::VectorXd::Zero(2);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
            es.compute(overlapmatrix);
            Eigen::MatrixXd sqrtm1 = es.operatorInverseSqrt();
            Eigen::MatrixXd prodmat = sqrtm1 * densitymatrix*sqrtm1;

            for (int i = 0; i < frag; i++) {
                fragmentCharges(0) += prodmat(i, i);
            }
            for (int i = frag; i < overlapmatrix.rows(); i++) {
                fragmentCharges(1) += prodmat(i, i);
            }

            return fragmentCharges;
        }

        std::vector<double> Orbitals::Oscillatorstrengths() const{
            std::vector<double> oscs;
            int size = _transition_dipoles.size();
            if (size > _BSE_singlet_energies.size()) {
                size = _BSE_singlet_energies.size();
            }
            for (int i = 0; i < size; ++i) {
                double osc = (_transition_dipoles[i] * _transition_dipoles[i]) * 2.0 / 3.0 * (_BSE_singlet_energies(i));
                oscs.push_back(osc);
            }
            return oscs;
        }

        double Orbitals::getTotalStateEnergy(const QMState& state)const{
          double total_energy=getQMEnergy()* tools::conv::ev2hrt;
          if (state.Type()==QMStateType::Gstate){
            return total_energy;
          }
          total_energy+=getExcitedStateEnergy(state);
          return total_energy;
        }
        
        double Orbitals::getExcitedStateEnergy(const QMState& state) const{

            double omega = 0.0;
            if(state.isTransition()){
              throw std::runtime_error("Total Energy does not exist for transition state");
            }

            if (state.Type() == QMStateType::Singlet) {
              if(BSESingletEnergies().size()<state.Index()+1){
                throw std::runtime_error("Orbitals::getTotalEnergy You want "+ state.ToString()+" which has not been calculated");
              }
                omega = BSESingletEnergies()[state.Index()];
            } else if (state.Type() == QMStateType::Triplet) {
               if(BSETripletEnergies().size()<state.Index()+1){
                  throw std::runtime_error("Orbitals::getTotalEnergy You want "+ state.ToString()+" which has not been calculated");
              }
                omega = BSETripletEnergies()[state.Index()];
            } else if (state.Type() == QMStateType::DQPstate) {
               if(this->QPdiagEnergies().size()<state.Index()+1-getGWAmin()){
                  throw std::runtime_error("Orbitals::getTotalEnergy You want "+ state.ToString()+" which has not been calculated");
              }
               return QPdiagEnergies()[state.Index()-getGWAmin()];
            }else if (state.Type() == QMStateType::KSstate) {
               if(this->MOEnergies().size()<state.Index()+1){
                  throw std::runtime_error("Orbitals::getTotalEnergy You want "+ state.ToString()+" which has not been calculated");
              }
               return QPdiagEnergies()[state.Index()];
            }else if (state.Type() == QMStateType::PQPstate) {
               if(this->_QPpert_energies.rows()<state.Index()+1-getGWAmin()){
                  throw std::runtime_error("Orbitals::getTotalEnergy You want "+ state.ToString()+" which has not been calculated");
              }
               return _QPpert_energies(state.Index()-getGWAmin(),3);
            } else {
                throw std::runtime_error("GetTotalEnergy only knows states:singlet,triplet,KS,DQP,PQP");
            }
            return  omega; //  e.g. hartree
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
        
        
        void Orbitals::OrderMOsbyEnergy(){
          std::vector<int> sort_index = this->SortEnergies();
          
          Eigen::MatrixXd MOcopy=this->MOCoefficients();
          Eigen::VectorXd Energy=this->MOEnergies();
          
          for(int i=0;i<Energy.size();++i){
            this->MOEnergies()(i)=Energy(sort_index[i]);
          }
          for(int i=0;i<Energy.size();++i){
            this->MOCoefficients().col(i)=MOcopy.col(sort_index[i]);
          }
        }

                /**
         * \brief Guess for a dimer based on monomer orbitals
         *
         * Given two monomer orbitals (A and B) constructs a guess for dimer
         * orbitals: | A 0 | and energies: [EA, EB]
         *           | 0 B |
         */
        void Orbitals::PrepareDimerGuess(const Orbitals& orbitalsA,const Orbitals& orbitalsB) {

            // constructing the direct product orbA x orbB
            int basisA = orbitalsA.getBasisSetSize();
            int basisB = orbitalsB.getBasisSetSize();

            int levelsA = orbitalsA.getNumberOfLevels();
            int levelsB = orbitalsB.getNumberOfLevels();

            int electronsA = orbitalsA.getNumberOfElectrons();
            int electronsB = orbitalsB.getNumberOfElectrons();

            
            this->MOCoefficients() = Eigen::MatrixXd::Zero(basisA + basisB, levelsA + levelsB);

            // AxB = | A 0 |  //   A = [EA, EB]  //
            //       | 0 B |  //                 //
            if(orbitalsA.getDFTbasis()!=orbitalsB.getDFTbasis()){
              throw std::runtime_error("Basissets of Orbitals A and B differ "+orbitalsA.getDFTbasis()+":"+orbitalsB.getDFTbasis());
            }
            this->setDFTbasis(orbitalsA.getDFTbasis());
            if(orbitalsA.getECP()!=orbitalsB.getECP()){
              throw std::runtime_error("ECPs of Orbitals A and B differ "+orbitalsA.getECP()+":"+orbitalsB.getECP());
            }
            this->setECP(orbitalsA.getECP());
            this->setBasisSetSize(basisA + basisB);
            this->setNumberOfLevels(electronsA + electronsB,
                    levelsA + levelsB - electronsA - electronsB);
            this->setNumberOfElectrons(electronsA + electronsB);

            this->MOCoefficients().block(0, 0, basisA, levelsA) = orbitalsA.MOCoefficients();
            this->MOCoefficients().block(basisA, levelsA, basisB, levelsB) = orbitalsB.MOCoefficients();

            Eigen::VectorXd& energies = this->MOEnergies();
            energies.resize(levelsA + levelsB);

            energies.segment(0, levelsA) = orbitalsA.MOEnergies();
            energies.segment(levelsA, levelsB) = orbitalsB.MOEnergies();
            
            this->OrderMOsbyEnergy();
            
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
                    if(split.size()<4){continue;}
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
            CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
            WriteToCpt(cpf);
        }

        void Orbitals::WriteToCpt(CheckpointFile f) const{
            WriteToCpt(f.getWriter("/QMdata"));
        }

        void Orbitals::WriteToCpt(CheckpointWriter w) const{
            w(XtpVersionStr(), "Version");
            w(_basis_set_size, "basis_set_size");
            w(_occupied_levels, "occupied_levels");
            w(_unoccupied_levels, "unoccupied_levels");
            w(_number_of_electrons, "number_of_electrons");

            w(_mo_energies, "mo_energies");
            w(_mo_coefficients, "mo_coefficients");

            // write qmatoms

            {
                CheckpointWriter qmasWriter = w.openChild("qmatoms");
                for (size_t idx = 0; idx < _atoms.size(); ++idx) {
                    auto qmaWriter =
                        qmasWriter.openChild("atom" + std::to_string(idx));
                    _atoms[idx]->WriteToCpt(qmaWriter);
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

            w(_useTDA, "useTDA");
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
        }

        void Orbitals::ReadFromCpt(const std::string& filename) {
            CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
            ReadFromCpt(cpf);
        }

        void Orbitals::ReadFromCpt(CheckpointFile f) {
            ReadFromCpt(f.getReader("/QMdata"));
        }

        void Orbitals::ReadFromCpt(CheckpointReader r) {
            r(_basis_set_size, "basis_set_size");
            r(_occupied_levels, "occupied_levels");
            r(_unoccupied_levels, "unoccupied_levels");
            r(_number_of_electrons, "number_of_electrons");

            r(_mo_energies, "mo_energies");
            r(_mo_coefficients, "mo_coefficients");

            // Read qmatoms
            {
                CheckpointReader qmasReader = r.openChild("qmatoms");
                size_t count = qmasReader.getNumDataSets();
                if(this->QMAtoms().size()>0){
                    std::vector< QMAtom* >::iterator it;
                    for (it = _atoms.begin(); it != _atoms.end(); ++it) delete *it;
                    _atoms.clear();
                }

                for (size_t i = 0; i < count; ++i) {
                    CheckpointReader qmaReader =
                        qmasReader.openChild("atom" + std::to_string(i));
                    QMAtom temp;
                    temp.ReadFromCpt(qmaReader);
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
            r(_useTDA, "useTDA");
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
        }
    }
}
