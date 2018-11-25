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

#ifndef _VOTCA_XTP_BSE_H
#define _VOTCA_XTP_BSE_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/qmstate.h>

namespace votca {
namespace xtp {

class BSE {

 public:
 
  BSE(Orbitals& orbitals,ctp::Logger &log,double min_print_weight):
        _log(log),
        _orbitals(orbitals),
        _eh_s(orbitals.eh_s()),
        _eh_t(orbitals.eh_t()),
        _bse_singlet_energies(orbitals.BSESingletEnergies()),
        _bse_singlet_coefficients(orbitals.BSESingletCoefficients()),
        _bse_singlet_coefficients_AR(orbitals.BSESingletCoefficientsAR()),
        _bse_triplet_energies(orbitals.BSETripletEnergies()),
        _bse_triplet_coefficients(orbitals.BSETripletCoefficients()),
        _min_print_weight(min_print_weight){};
  
  void setGWData(const TCMatrix_gwbse* Mmn,const PPM* ppm,const Eigen::MatrixXd* Hqp){
      _Mmn=Mmn;
      _ppm=ppm;
      _Hqp=Hqp;   
  }
  
  void setBSEindices(int homo,int vmin, int cmax, int nmax) {
                _homo=homo;
                _bse_vmin = vmin;
                _bse_vmax = homo;
                _bse_cmin = homo+1;
                _bse_cmax = cmax;
                _bse_nmax = nmax;
                _bse_vtotal = _bse_vmax - _bse_vmin + 1;
                _bse_ctotal = _bse_cmax - _bse_cmin + 1;
                _bse_size = _bse_vtotal * _bse_ctotal;
                // indexing info BSE vector index to occupied/virtual orbital
                for (int v = 0; v < _bse_vtotal; v++) {
                    for (int c = 0; c < _bse_ctotal; c++) {
                        _index2v.push_back(_bse_vmin + v);
                        _index2c.push_back(_bse_cmin + c);
                    }
                }
                return;
            }

   
  void Solve_triplets();
  void Solve_singlets();
  void Solve_singlets_BTDA();
  void Analyze_triplets(const AOBasis& dftbasis);
  void Analyze_singlets(const AOBasis& dftbasis);
   
  void FreeMatrices(){
      _eh_t.resize(0, 0);
      _eh_s.resize(0, 0);
  }
  
  void SetupHs();
  
  void SetupHt();
  
  void FreeTriplets(){
      _bse_triplet_coefficients.resize(0,0);
  }
  
  void FreeSinglets(){
      _bse_singlet_coefficients.resize(0,0);
      _bse_singlet_coefficients_AR.resize(0,0);
  }
 
 private:

     struct Interaction {
    Eigen::VectorXd exchange_contrib;
    Eigen::VectorXd direct_contrib;
    Eigen::VectorXd qp_contrib;
};

struct Population {
    std::vector<Eigen::VectorXd> popH;
    std::vector<Eigen::VectorXd> popE;
    std::vector<Eigen::VectorXd> Crgs;
    Eigen::VectorXd popGs;
};   
 
      
ctp::Logger &_log;
  int  _homo;
  int  _bse_vmin;
  int  _bse_vmax;
  int  _bse_cmin;
  int  _bse_cmax;
  int  _bse_size;
  int  _bse_vtotal;
  int  _bse_ctotal;
  int _bse_nmax;
  
  Orbitals& _orbitals;
  
  const TCMatrix_gwbse* _Mmn;
  const PPM* _ppm;
  const Eigen::MatrixXd* _Hqp;
  

  // BSE variables and functions
  MatrixXfd& _eh_s;  // only for storage in orbitals object
  MatrixXfd& _eh_t;  // only for storage in orbitals object

  VectorXfd& _bse_singlet_energies;  // stored in orbitals object
  MatrixXfd& _bse_singlet_coefficients;  // stored in orbitals
                                                      // object
  MatrixXfd& _bse_singlet_coefficients_AR;  // stored in orbitals
                                                         // object
  VectorXfd& _bse_triplet_energies;  // stored in orbitals object
  MatrixXfd& _bse_triplet_coefficients;  // stored in orbitals
                                                      // object
  
  double _min_print_weight;

   template <typename T>
  void Add_Hqp(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H);
   template <typename T>
  void Add_Hx(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor);
   template <typename T>
   void Add_Hd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H);
   template <typename T>
  void Add_Hd2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor);
  
  std::vector<int> _index2v;
  std::vector<int> _index2c;

 void printFragInfo(const Population& pop, int i);
 void printWeights(int i_bse, double weight);
 
  Interaction Analyze_eh_interaction(const QMStateType& type);
  Eigen::VectorXd Analyze_IndividualContribution(const QMStateType& type, const MatrixXfd& H);

  Population FragmentPopulations(const QMStateType& type, const AOBasis& dftbasis);

  std::vector<Eigen::MatrixXd > CalcFreeTransition_Dipoles(const AOBasis& dftbasis);

  std::vector<tools::vec > CalcCoupledTransition_Dipoles(const AOBasis& dftbasis);
};
}
}

#endif /* _VOTCA_XTP_BSE_H */
