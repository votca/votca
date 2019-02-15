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

#ifndef __XTP_NUMERICAL_INTEGRATION__H
#define __XTP_NUMERICAL_INTEGRATION__H

#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/grid_containers.h>
#include <votca/xtp/gridbox.h>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/vxc_functionals.h>

#include <xc.h>
#undef LOG

namespace votca {
namespace xtp {
class LebedevGrid;

struct Gyrationtensor {
  double mass;
  tools::vec centroid;
  tools::matrix gyration;
};

class NumericalIntegration {
 public:
  NumericalIntegration() : _density_set(false), _setXC(false){};

  ~NumericalIntegration();

  void GridSetup(const std::string& type, std::vector<QMAtom*> atoms,
                 const AOBasis& basis);

  double getExactExchange(const std::string& functional);
  std::vector<const tools::vec*> getGridpoints() const;
  std::vector<double> getWeightedDensities() const;
  int getGridSize() const { return _totalgridsize; }
  unsigned getBoxesSize() const { return _grid_boxes.size(); }

  void setXCfunctional(const std::string& functional);
  double IntegrateDensity(const Eigen::MatrixXd& density_matrix);
  double IntegratePotential(const tools::vec& rvector);
  Eigen::MatrixXd IntegratePotential(const AOBasis& externalbasis);

  Eigen::MatrixXd IntegrateExternalPotential(
      const std::vector<double>& Potentialvalues);
  Gyrationtensor IntegrateGyrationTensor(const Eigen::MatrixXd& density_matrix);
  Eigen::MatrixXd IntegrateVXC(const Eigen::MatrixXd& density_matrix);
  double getTotEcontribution() { return _EXC; }

 private:
  void FindSignificantShells(const AOBasis& basis);
  void EvaluateXC(const double rho, const double sigma, double& f_xc,
                  double& df_drho, double& df_dsigma);
  double erf1c(double x);

  void SortGridpointsintoBlocks(
      std::vector<std::vector<GridContainers::Cartesian_gridpoint> >& grid);

  Eigen::MatrixXd CalcInverseAtomDist(std::vector<QMAtom*>& atoms);
  int UpdateOrder(LebedevGrid& sphericalgridofElement, int maxorder,
                  std::vector<double>& PruningIntervals, double r);

  GridContainers::Cartesian_gridpoint CreateCartesianGridpoint(
      const tools::vec& atomA_pos, GridContainers::radial_grid& radial_grid,
      GridContainers::spherical_grid& spherical_grid, unsigned i_rad,
      unsigned i_sph);

  Eigen::VectorXd SSWpartition(int igrid, const Eigen::MatrixXd& rq,
                               const Eigen::MatrixXd& Rij);
  void SSWpartitionAtom(
      std::vector<QMAtom*>& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid,
      unsigned i_atom, const Eigen::MatrixXd& Rij);
  Eigen::MatrixXd CalcDistanceAtomsGridpoints(
      std::vector<QMAtom*>& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid);

  int _totalgridsize;
  std::vector<GridBox> _grid_boxes;
  std::vector<unsigned> thread_start;
  std::vector<unsigned> thread_stop;
  int xfunc_id;
  double _EXC;
  bool _density_set;
  bool _setXC;
  int _AOBasisSize;

  bool _use_separate;
  int cfunc_id;
  xc_func_type xfunc;  // handle for exchange functional
  xc_func_type cfunc;  // handle for correlation functional
};

}  // namespace xtp
}  // namespace votca
#endif /* NUMERICAL_INTEGRATION_H */
