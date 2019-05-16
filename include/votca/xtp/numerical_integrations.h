/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef XTP_NUMERICAL_INTEGRATION_H
#define XTP_NUMERICAL_INTEGRATION_H

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/grid_containers.h>
#include <votca/xtp/gridbox.h>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/vxc_functionals.h>

#include <xc.h>
#undef LOG

namespace votca {
namespace xtp {
class LebedevGrid;

struct Gyrationtensor {
  double mass;
  Eigen::Vector3d centroid;
  Eigen::Matrix3d gyration;
};

class NumericalIntegration {
 public:
  ~NumericalIntegration();

  void GridSetup(const std::string& type, const QMMolecule& atoms,
                 const AOBasis& basis);

  double getExactExchange(const std::string& functional) const;
  std::vector<const Eigen::Vector3d*> getGridpoints() const;
  std::vector<double> getWeightedDensities() const;
  int getGridSize() const { return _totalgridsize; }
  int getBoxesSize() const { return _grid_boxes.size(); }

  void setXCfunctional(const std::string& functional);
  double IntegrateDensity(const Eigen::MatrixXd& density_matrix);
  double IntegratePotential(const Eigen::Vector3d& rvector) const;
  Eigen::Vector3d IntegrateField(const Eigen::Vector3d& rvector) const;
  Eigen::MatrixXd IntegratePotential(const AOBasis& externalbasis) const;

  Gyrationtensor IntegrateGyrationTensor(const Eigen::MatrixXd& density_matrix);

  struct E_Vxc {
    Eigen::MatrixXd Vxc;
    double Exc;
  };

  Mat_p_Energy IntegrateVXC(const Eigen::MatrixXd& density_matrix);

 private:
  struct XC_entry {
    double f_xc = 0;  // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r
    double df_drho = 0;    // v_xc_rho(r) = df/drho
    double df_dsigma = 0;  // df/dsigma ( df/dgrad(rho) = df/dsigma *
                           // dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
  };

  Eigen::VectorXd CalcAOValue_and_Grad(Eigen::MatrixX3d& ao_grad,
                                       const GridBox& box,
                                       const Eigen::Vector3d& point) const;
  Eigen::VectorXd CalcAOValues(const GridBox& box,
                               const Eigen::Vector3d& pos) const;
  void FindSignificantShells(const AOBasis& basis);
  XC_entry EvaluateXC(double rho, double sigma) const;
  double erf1c(double x) const;

  void SortGridpointsintoBlocks(
      std::vector<std::vector<GridContainers::Cartesian_gridpoint> >& grid);

  Eigen::MatrixXd CalcInverseAtomDist(const QMMolecule& atoms) const;
  int UpdateOrder(LebedevGrid& sphericalgridofElement, int maxorder,
                  std::vector<double>& PruningIntervals, double r) const;

  GridContainers::Cartesian_gridpoint CreateCartesianGridpoint(
      const Eigen::Vector3d& atomA_pos,
      GridContainers::radial_grid& radial_grid,
      GridContainers::spherical_grid& spherical_grid, int i_rad,
      int i_sph) const;

  Eigen::VectorXd SSWpartition(const Eigen::VectorXd& rq_i,
                               const Eigen::MatrixXd& Rij) const;
  void SSWpartitionAtom(
      const QMMolecule& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid, int i_atom,
      const Eigen::MatrixXd& Rij) const;
  Eigen::MatrixXd CalcDistanceAtomsGridpoints(
      const QMMolecule& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid) const;

  int _totalgridsize;
  std::vector<GridBox> _grid_boxes;
  std::vector<unsigned> thread_start;
  std::vector<unsigned> thread_stop;
  int xfunc_id;
  bool _density_set = false;
  bool _setXC = false;

  bool _use_separate;
  int cfunc_id;
  xc_func_type xfunc;  // handle for exchange functional
  xc_func_type cfunc;  // handle for correlation functional
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_NUMERICAL_INTEGRATION_H
