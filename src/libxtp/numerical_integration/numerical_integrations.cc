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

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iterator>
#include <string>
#include <votca/tools/constants.h>
#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/vxc_functionals.h>

namespace votca {
namespace xtp {

NumericalIntegration::~NumericalIntegration() {
  if (_setXC) {
    xc_func_end(&xfunc);
    if (_use_separate) {
      xc_func_end(&cfunc);
    }
  }
}

double NumericalIntegration::getExactExchange(
    const std::string& functional) const {

  double exactexchange = 0.0;
  Vxc_Functionals map;
  std::vector<std::string> strs;

  boost::split(strs, functional, boost::is_any_of(" "));
  if (strs.size() > 2) {
    throw std::runtime_error("Too many functional names");
  } else if (strs.size() < 1) {
    throw std::runtime_error("Specify at least one functional");
  }

  for (unsigned i = 0; i < strs.size(); i++) {

    int func_id = map.getID(strs[i]);
    if (func_id < 0) {
      exactexchange = 0.0;
      break;
    }
    xc_func_type func;
    if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
      throw std::runtime_error(
          (boost::format("Functional %s not found\n") % strs[i]).str());
    }
    if (exactexchange > 0 && func.cam_alpha > 0) {
      throw std::runtime_error(
          "You have specified two functionals with exact exchange");
    }
    exactexchange += func.cam_alpha;
    xc_func_end(&func);
  }

  return exactexchange;
}

void NumericalIntegration::setXCfunctional(const std::string& functional) {

  Vxc_Functionals map;
  std::vector<std::string> strs;
  tools::Tokenizer tok(functional, " ,\n\t");
  tok.ToVector(strs);
  xfunc_id = 0;
  _use_separate = false;
  cfunc_id = 0;
  if (strs.size() == 1) {
    xfunc_id = map.getID(strs[0]);
  } else if (strs.size() == 2) {
    xfunc_id = map.getID(strs[0]);
    cfunc_id = map.getID(strs[1]);
    _use_separate = true;
  } else {
    std::cout << "LIBXC " << strs.size() << std::endl;
    throw std::runtime_error(
        "LIBXC. Please specify one combined or an exchange and a correlation "
        "functionals");
  }

  if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
    throw std::runtime_error(
        (boost::format("Functional %s not found\n") % strs[0]).str());
  }
  if (xfunc.info->kind != 2 && !_use_separate) {
    throw std::runtime_error(
        "Your functional misses either correlation or exchange, please specify "
        "another functional, separated by whitespace");
  }
  if (_use_separate) {
    if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
      throw std::runtime_error(
          (boost::format("Functional %s not found\n") % strs[1]).str());
    }
    if ((xfunc.info->kind + cfunc.info->kind) != 1) {
      throw std::runtime_error(
          "Your functionals are not one exchange and one correlation");
    }
  }
  _setXC = true;
  return;
}

NumericalIntegration::XC_entry NumericalIntegration::EvaluateXC(
    double rho, double sigma) const {

  NumericalIntegration::XC_entry result;
  switch (xfunc.info->family) {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&xfunc, 1, &rho, &result.f_xc, &result.df_drho);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, &result.f_xc, &result.df_drho,
                     &result.df_dsigma);
      break;
  }
  if (_use_separate) {
    NumericalIntegration::XC_entry temp;
    // via libxc correlation part only
    switch (cfunc.info->family) {
      case XC_FAMILY_LDA:
        xc_lda_exc_vxc(&cfunc, 1, &rho, &temp.f_xc, &temp.df_drho);
        break;
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA:
        xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, &temp.f_xc, &temp.df_drho,
                       &temp.df_dsigma);
        break;
    }

    result.f_xc += temp.f_xc;
    result.df_drho += temp.df_drho;
    result.df_dsigma += temp.df_dsigma;
  }

  return result;
}

double NumericalIntegration::IntegratePotential(
    const Eigen::Vector3d& rvector) const {

  double result = 0.0;
  assert(_density_set && "Density not calculated");
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid_boxes[i].getGridPoints();
    const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
    const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
    for (unsigned j = 0; j < points.size(); j++) {
      double charge = -weights[j] * densities[j];
      double dist = (points[j] - rvector).norm();
      result = charge / dist;
    }
  }
  return result;
}

Vector9d NumericalIntegration::IntegrateV(
    const Eigen::Vector3d& rvector) const {

  Vector9d result = Vector9d::Zero();
  assert(_density_set && "Density not calculated");
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid_boxes[i].getGridPoints();
    const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
    const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
    for (unsigned j = 0; j < points.size(); j++) {
      double charge = -weights[j] * densities[j];
      Eigen::Vector3d r = points[j] - rvector;
      double dist = r.norm();
      result(0) += charge / dist;
      double dist2 = dist * dist;
      double dist3 = dist * dist2;
      result.segment<3>(1) += charge * r / dist3;  // x,y,z
      double charge_dist5 = charge / std::pow(dist, 5);
      result(4) += charge_dist5 / 2 * (9 * r.z() * r.z() - 3 * dist2);
      double fac = 3 * std::sqrt(3) * charge_dist5;
      result(5) += fac * r.x() * r.z();
      result(6) += fac * r.y() * r.z();
      result(7) += fac / 2 * (r.x() * r.x() - r.y() * r.y());
      result(8) += fac * r.x() * r.y();
    }
  }
  return result;
}

void NumericalIntegration::SortGridpointsintoBlocks(
    std::vector<std::vector<GridContainers::Cartesian_gridpoint> >& grid) {
  const double boxsize = 1;  // 1 bohr

  Eigen::Vector3d min =
      Eigen::Vector3d::Ones() * std::numeric_limits<double>::max();
  Eigen::Vector3d max =
      Eigen::Vector3d::Ones() * std::numeric_limits<double>::min();

  for (unsigned i = 0; i < grid.size(); i++) {
    for (unsigned j = 0; j < grid[i].size(); j++) {
      const Eigen::Vector3d& pos = grid[i][j].grid_pos;
      if (pos[0] > max[0]) {
        max[0] = pos[0];
      } else if (pos[0] < min[0]) {
        min[0] = pos[0];
      }
      if (pos[1] > max[1]) {
        max[1] = pos[1];
      } else if (pos[1] < min[1]) {
        min[1] = pos[1];
      }
      if (pos[2] > max[2]) {
        max[2] = pos[2];
      } else if (pos[2] < min[2]) {
        min[2] = pos[2];
      }
    }
  }

  Eigen::Vector3d molextension = max - min;
  Eigen::Vector3d numberofboxes = molextension / boxsize;
  Eigen::Vector3d roundednumofbox(std::ceil(numberofboxes[0]),
                                  std::ceil(numberofboxes[1]),
                                  std::ceil(numberofboxes[2]));

  std::vector<std::vector<
      std::vector<std::vector<GridContainers::Cartesian_gridpoint*> > > >
      boxes;
  // creating temparray
  for (unsigned i = 0; i < unsigned(roundednumofbox[0]); i++) {
    std::vector<
        std::vector<std::vector<GridContainers::Cartesian_gridpoint*> > >
        boxes_yz;
    for (unsigned j = 0; j < unsigned(roundednumofbox[1]); j++) {
      std::vector<std::vector<GridContainers::Cartesian_gridpoint*> > boxes_z;
      for (unsigned k = 0; k < unsigned(roundednumofbox[2]); k++) {
        std::vector<GridContainers::Cartesian_gridpoint*> box;
        box.reserve(100);
        boxes_z.push_back(box);
      }
      boxes_yz.push_back(boxes_z);
    }
    boxes.push_back(boxes_yz);
  }

  for (auto& atomgrid : grid) {
    for (auto& gridpoint : atomgrid) {
      Eigen::Vector3d pos = gridpoint.grid_pos - min;
      Eigen::Vector3d index = pos / boxsize;
      int i_x = int(index[0]);
      int i_y = int(index[1]);
      int i_z = int(index[2]);
      boxes[i_x][i_y][i_z].push_back(&gridpoint);
    }
  }

  for (auto& boxes_xy : boxes) {
    for (auto& boxes_z : boxes_xy) {
      for (auto& box : boxes_z) {
        if (box.size() < 1) {
          continue;
        }
        GridBox gridbox;

        for (const auto& point : box) {
          gridbox.addGridPoint(*point);
        }
        _grid_boxes.push_back(gridbox);
      }
    }
  }
  return;
}

void NumericalIntegration::FindSignificantShells(const AOBasis& basis) {
  for (unsigned i = 0; i < _grid_boxes.size(); ++i) {
    GridBox& box = _grid_boxes[i];
    for (const AOShell& store : basis) {
      const double decay = store.getMinDecay();
      const Eigen::Vector3d& shellpos = store.getPos();
      for (const auto& point : box.getGridPoints()) {
        Eigen::Vector3d dist = shellpos - point;
        double distsq = dist.squaredNorm();
        // if contribution is smaller than -ln(1e-10), add shell to list
        if ((decay * distsq) < 20.7) {
          box.addShell(&store);
          break;
        }
      }
    }
  }

  std::vector<GridBox> grid_boxes_copy;
  int combined = 0;
  // use vector of bool to indicate if a gridbox has already been merged into
  // another
  std::vector<bool> Merged = std::vector<bool>(_grid_boxes.size(), false);
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    if (Merged[i]) {
      continue;
    }
    GridBox box = _grid_boxes[i];
    if (box.Shellsize() < 1) {
      continue;
    }
    Merged[i] = true;
    for (unsigned j = i + 1; j < _grid_boxes.size(); j++) {
      if (GridBox::compareGridboxes(_grid_boxes[i], _grid_boxes[j])) {
        Merged[j] = true;
        box.addGridBox(_grid_boxes[j]);
        combined++;
      }
    }
    grid_boxes_copy.push_back(box);
  }

  _grid_boxes = grid_boxes_copy;
  for (auto& box : _grid_boxes) {
    box.PrepareForIntegration();
  }
}

Eigen::VectorXd NumericalIntegration::CalcAOValue_and_Grad(
    Eigen::MatrixX3d& ao_grad, const GridBox& box,
    const Eigen::Vector3d& point) const {
  Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
  const std::vector<GridboxRange>& aoranges = box.getAOranges();
  const std::vector<const AOShell*>& shells = box.getShells();
  for (unsigned j = 0; j < box.Shellsize(); ++j) {
    Eigen::Block<Eigen::MatrixX3d> grad_block =
        ao_grad.block(aoranges[j].start, 0, aoranges[j].size, 3);
    Eigen::VectorBlock<Eigen::VectorXd> ao_block =
        ao.segment(aoranges[j].start, aoranges[j].size);
    shells[j]->EvalAOspace(ao_block, grad_block, point);
  }
  return ao;
}

Mat_p_Energy NumericalIntegration::IntegrateVXC(
    const Eigen::MatrixXd& density_matrix) {

  int nthreads = OPENMP::getMaxThreads();

  std::vector<Eigen::MatrixXd> vxc_thread = std::vector<Eigen::MatrixXd>(
      nthreads,
      Eigen::MatrixXd::Zero(density_matrix.rows(), density_matrix.cols()));
  std::vector<double> Exc_thread = std::vector<double>(nthreads, 0.0);

#pragma omp parallel for schedule(guided)
  for (unsigned i = 0; i < _grid_boxes.size(); ++i) {

    double EXC_box = 0.0;
    const GridBox& box = _grid_boxes[i];
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const Eigen::MatrixXd DMAT_symm = DMAT_here + DMAT_here.transpose();
    double cutoff = 1.e-40 / density_matrix.rows() / density_matrix.rows();
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }
    Eigen::MatrixXd Vxc_here =
        Eigen::MatrixXd::Zero(DMAT_here.rows(), DMAT_here.cols());
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    // iterate over gridpoints
    for (unsigned p = 0; p < box.size(); p++) {

      Eigen::MatrixX3d ao_grad = Eigen::MatrixX3d::Zero(box.Matrixsize(), 3);
      Eigen::VectorXd ao = CalcAOValue_and_Grad(ao_grad, box, points[p]);
      const double rho = 0.5 * (ao.transpose() * DMAT_symm * ao).value();
      const double weight = weights[p];
      if (rho * weight < 1.e-20)
        continue;  // skip the rest, if density is very small
      const Eigen::Vector3d rho_grad = ao.transpose() * DMAT_symm * ao_grad;
      const double sigma = (rho_grad.transpose() * rho_grad).value();
      const Eigen::VectorXd grad = ao_grad * rho_grad;
      NumericalIntegration::XC_entry xc = EvaluateXC(rho, sigma);
      EXC_box += weight * rho * xc.f_xc;
      auto addXC = weight * (0.5 * xc.df_drho * ao + 2.0 * xc.df_dsigma * grad);
      Vxc_here.noalias() += addXC * ao.transpose();
    }
    box.AddtoBigMatrix(vxc_thread[OPENMP::getThreadId()], Vxc_here);
    Exc_thread[OPENMP::getThreadId()] += EXC_box;
  }

  double EXC = std::accumulate(Exc_thread.begin(), Exc_thread.end(), 0.0);
  Eigen::MatrixXd Vxc = std::accumulate(
      vxc_thread.begin(), vxc_thread.end(),
      Eigen::MatrixXd::Zero(density_matrix.rows(), density_matrix.cols())
          .eval());

  Mat_p_Energy Oxc(EXC, Vxc + Vxc.transpose());
  return Oxc;
}

double NumericalIntegration::IntegrateDensity(
    const Eigen::MatrixXd& density_matrix) {

  int nthreads = OPENMP::getMaxThreads();
  std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);

#pragma omp parallel for schedule(guided)
  for (unsigned i = 0; i < _grid_boxes.size(); ++i) {
    double N_box = 0.0;
    GridBox& box = _grid_boxes[i];
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    box.prepareDensity();
    // iterate over gridpoints
    for (unsigned p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = CalcAOValues(box, points[p]);
      double rho = (ao.transpose() * DMAT_here * ao)(0, 0);
      box.addDensity(rho);
      N_box += rho * weights[p];
    }
    N_thread[OPENMP::getThreadId()] += N_box;
  }
  double N = std::accumulate(N_thread.begin(), N_thread.end(), 0.0);
  _density_set = true;
  return N;
}

Eigen::VectorXd NumericalIntegration::CalcAOValues(
    const GridBox& box, const Eigen::Vector3d& pos) const {
  const std::vector<GridboxRange>& aoranges = box.getAOranges();
  const std::vector<const AOShell*> shells = box.getShells();
  Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
  for (unsigned j = 0; j < box.Shellsize(); ++j) {
    Eigen::VectorBlock<Eigen::VectorXd> ao_block =
        ao.segment(aoranges[j].start, aoranges[j].size);
    shells[j]->EvalAOspace(ao_block, pos);
  }
  return ao;
}

Gyrationtensor NumericalIntegration::IntegrateGyrationTensor(
    const Eigen::MatrixXd& density_matrix) {

  int nthreads = OPENMP::getMaxThreads();
  std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);
  std::vector<Eigen::Vector3d> centroid_thread =
      std::vector<Eigen::Vector3d>(nthreads, Eigen::Vector3d::Zero());
  std::vector<Eigen::Matrix3d> gyration_thread =
      std::vector<Eigen::Matrix3d>(nthreads, Eigen::Matrix3d::Zero());

#pragma omp parallel for schedule(guided)
  for (unsigned i = 0; i < _grid_boxes.size(); ++i) {
    double N_box = 0.0;
    Eigen::Vector3d centroid_box = Eigen::Vector3d::Zero();
    Eigen::Matrix3d gyration_box = Eigen::Matrix3d::Zero();
    GridBox& box = _grid_boxes[i];
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    box.prepareDensity();
    // iterate over gridpoints
    for (unsigned p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = CalcAOValues(box, points[p]);
      double rho = (ao.transpose() * DMAT_here * ao).value();
      box.addDensity(rho);
      N_box += rho * weights[p];
      centroid_box += rho * weights[p] * points[p];
      gyration_box += rho * weights[p] * points[p] * points[p].transpose();
    }
    N_thread[OPENMP::getThreadId()] += N_box;
    centroid_thread[OPENMP::getThreadId()] += centroid_box;
    gyration_thread[OPENMP::getThreadId()] += gyration_box;
  }
  double N = std::accumulate(N_thread.begin(), N_thread.end(), 0.0);
  Eigen::Vector3d centroid =
      std::accumulate(centroid_thread.begin(), centroid_thread.end(),
                      Eigen::Vector3d::Zero().eval());
  Eigen::Matrix3d gyration =
      std::accumulate(gyration_thread.begin(), gyration_thread.end(),
                      Eigen::Matrix3d::Zero().eval());

  _density_set = true;
  // Normalize
  centroid = centroid / N;
  gyration = gyration / N;
  gyration = gyration - centroid * centroid.transpose();
  Gyrationtensor gyro;
  gyro.mass = N;
  gyro.centroid = centroid;
  gyro.gyration = gyration;

  return gyro;
}

std::vector<const Eigen::Vector3d*> NumericalIntegration::getGridpoints()
    const {
  std::vector<const Eigen::Vector3d*> gridpoints;
  gridpoints.reserve(this->getGridSize());
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid_boxes[i].getGridPoints();
    for (unsigned j = 0; j < points.size(); j++) {
      gridpoints.push_back(&points[j]);
    }
  }
  return gridpoints;
}

Eigen::MatrixXd NumericalIntegration::IntegratePotential(
    const AOBasis& externalbasis) const {
  Eigen::MatrixXd Potential = Eigen::MatrixXd::Zero(
      externalbasis.AOBasisSize(), externalbasis.AOBasisSize());

  assert(_density_set && "Density not calculated");
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid_boxes[i].getGridPoints();
    const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
    const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
    for (unsigned j = 0; j < points.size(); j++) {
      double weighteddensity = weights[j] * densities[j];
      if (weighteddensity < 1e-12) {
        continue;
      }
      AOESP esp;
      esp.setPosition(points[j]);
      esp.Fill(externalbasis);
      Potential += weighteddensity * esp.Matrix();
    }
  }
  return Potential;
}

Eigen::MatrixXd NumericalIntegration::CalcInverseAtomDist(
    const QMMolecule& atoms) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(atoms.size(), atoms.size());
#pragma omp parallel for
  for (int i = 0; i < atoms.size(); ++i) {
    const Eigen::Vector3d& pos_a = atoms[i].getPos();
    for (int j = 0; j < i; ++j) {
      const Eigen::Vector3d& pos_b = atoms[j].getPos();
      result(j, i) = 1 / (pos_a - pos_b).norm();
    }
  }
  return result + result.transpose();
}

int NumericalIntegration::UpdateOrder(LebedevGrid& sphericalgridofElement,
                                      int maxorder,
                                      std::vector<double>& PruningIntervals,
                                      double r) const {
  int order;
  int maxindex = sphericalgridofElement.getIndexFromOrder(maxorder);
  if (maxindex == 1) {
    // smallest possible grid anyway, nothing to do
    order = maxorder;
  } else if (maxindex == 2) {
    // only three intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(1);  // 1;
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[3])) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else {
      order = sphericalgridofElement.getOrderFromIndex(1);
    }  // maxorder == 2
  } else {
    // five intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[1])) {
      order = sphericalgridofElement.getOrderFromIndex(4);
    } else if ((r >= PruningIntervals[1]) && (r < PruningIntervals[2])) {
      order =
          sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 4));
    } else if ((r >= PruningIntervals[2]) && (r < PruningIntervals[3])) {
      order = maxorder;
    } else {
      order =
          sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 1));
    }
  }
  return order;
}

GridContainers::Cartesian_gridpoint
    NumericalIntegration::CreateCartesianGridpoint(
        const Eigen::Vector3d& atomA_pos,
        GridContainers::radial_grid& radial_grid,
        GridContainers::spherical_grid& spherical_grid, int i_rad,
        int i_sph) const {
  GridContainers::Cartesian_gridpoint gridpoint;
  double p = spherical_grid.phi[i_sph];
  double t = spherical_grid.theta[i_sph];
  const Eigen::Vector3d s =
      Eigen::Vector3d{sin(p) * cos(t), sin(p) * sin(t), cos(p)};
  double r = radial_grid.radius[i_rad];
  gridpoint.grid_pos = atomA_pos + r * s;
  gridpoint.grid_weight =
      radial_grid.weight[i_rad] * spherical_grid.weight[i_sph];
  return gridpoint;
}

Eigen::MatrixXd NumericalIntegration::CalcDistanceAtomsGridpoints(
    const QMMolecule& atoms,
    std::vector<GridContainers::Cartesian_gridpoint>& atomgrid) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(atoms.size(), atomgrid.size());
#pragma omp parallel for
  for (int i = 0; i < atoms.size(); ++i) {
    const Eigen::Vector3d& atom_pos = atoms[i].getPos();
    for (unsigned j = 0; j < atomgrid.size(); ++j) {
      const auto& gridpoint = atomgrid[j];
      result(i, j) = (atom_pos - gridpoint.grid_pos).norm();
    }
  }
  return result;
}

void NumericalIntegration::SSWpartitionAtom(
    const QMMolecule& atoms,
    std::vector<GridContainers::Cartesian_gridpoint>& atomgrid, int i_atom,
    const Eigen::MatrixXd& Rij) const {
  Eigen::MatrixXd AtomGridDist = CalcDistanceAtomsGridpoints(atoms, atomgrid);

#pragma omp parallel for schedule(guided)
  for (int i_grid = 0; i_grid < int(atomgrid.size()); i_grid++) {
    Eigen::VectorXd p = SSWpartition(AtomGridDist.col(i_grid), Rij);
    // check weight sum
    double wsum = p.sum();
    if (wsum != 0.0) {
      // update the weight of this grid point
      atomgrid[i_grid].grid_weight *= p[i_atom] / wsum;
    } else {
      std::cerr << "\nSum of partition weights of grid point " << i_grid
                << " of atom " << i_atom << " is zero! ";
      throw std::runtime_error("\nThis should never happen!");
    }
  }  // partition weight for each gridpoint
}

void NumericalIntegration::GridSetup(const std::string& type,
                                     const QMMolecule& atoms,
                                     const AOBasis& basis) {
  GridContainers initialgrids;
  // get radial grid per element
  EulerMaclaurinGrid radialgridofElement;
  initialgrids.radial_grids = radialgridofElement.CalculateAtomicRadialGrids(
      basis, atoms, type);  // this checks out 1:1 with NWChem results! AWESOME
  LebedevGrid sphericalgridofElement;
  initialgrids.spherical_grids =
      sphericalgridofElement.CalculateSphericalGrids(atoms, type);

  // for the partitioning, we need all inter-center distances later, stored in
  // matrix
  Eigen::MatrixXd Rij = CalcInverseAtomDist(atoms);
  _totalgridsize = 0;
  std::vector<std::vector<GridContainers::Cartesian_gridpoint> > grid;

  for (int i_atom = 0; i_atom < atoms.size(); ++i_atom) {
    const QMAtom& atom = atoms[i_atom];

    const Eigen::Vector3d& atomA_pos = atom.getPos();
    const std::string& name = atom.getElement();
    GridContainers::radial_grid radial_grid =
        initialgrids.radial_grids.at(name);
    GridContainers::spherical_grid spherical_grid =
        initialgrids.spherical_grids.at(name);

    // maximum order (= number of points) in spherical integration grid
    int maxorder = sphericalgridofElement.Type2MaxOrder(name, type);
    // for pruning of integration grid, get interval boundaries for this element
    std::vector<double> PruningIntervals =
        radialgridofElement.CalculatePruningIntervals(name);
    int current_order = 0;
    // for each radial value
    std::vector<GridContainers::Cartesian_gridpoint> atomgrid;
    for (int i_rad = 0; i_rad < radial_grid.radius.size(); i_rad++) {
      double r = radial_grid.radius[i_rad];

      // which Lebedev order for this point?
      int order =
          UpdateOrder(sphericalgridofElement, maxorder, PruningIntervals, r);
      // get new spherical grid, if order changed
      if (order != current_order) {
        spherical_grid = sphericalgridofElement.CalculateUnitSphereGrid(order);
        current_order = order;
      }

      for (int i_sph = 0; i_sph < spherical_grid.phi.size(); i_sph++) {
        GridContainers::Cartesian_gridpoint gridpoint =
            CreateCartesianGridpoint(atomA_pos, radial_grid, spherical_grid,
                                     i_rad, i_sph);
        atomgrid.push_back(gridpoint);
      }  // spherical gridpoints
    }    // radial gridpoint

    SSWpartitionAtom(atoms, atomgrid, i_atom, Rij);
    // now remove points from the grid with negligible weights
    std::vector<GridContainers::Cartesian_gridpoint> atomgrid_cleanedup;
    for (const auto& point : atomgrid) {
      if (point.grid_weight > 1e-13) {
        atomgrid_cleanedup.push_back(point);
      }
    }

    _totalgridsize += atomgrid_cleanedup.size();
    grid.push_back(atomgrid_cleanedup);
  }  // atoms
  SortGridpointsintoBlocks(grid);
  FindSignificantShells(basis);
  return;
}

Eigen::VectorXd NumericalIntegration::SSWpartition(
    const Eigen::VectorXd& rq_i, const Eigen::MatrixXd& Rij) const {
  const double ass = 0.725;
  // initialize partition vector to 1.0
  Eigen::VectorXd p = Eigen::VectorXd::Ones(rq_i.size());
  const double tol_scr = 1e-10;
  const double leps = 1e-6;
  // go through centers
  for (int i = 1; i < rq_i.size(); i++) {
    double rag = rq_i(i);
    // through all other centers (one-directional)
    for (int j = 0; j < i; j++) {
      if ((std::abs(p[i]) > tol_scr) || (std::abs(p[j]) > tol_scr)) {
        double mu = (rag - rq_i(j)) * Rij(j, i);
        if (mu > ass) {
          p[i] = 0.0;
        } else if (mu < -ass) {
          p[j] = 0.0;
        } else {
          double sk;
          if (std::abs(mu) < leps) {
            sk = -1.88603178008 * mu + 0.5;
          } else {
            sk = erf1c(mu);
          }
          if (mu > 0.0) sk = 1.0 - sk;
          p[j] = p[j] * sk;
          p[i] = p[i] * (1.0 - sk);
        }
      }
    }
  }
  return p;
}

double NumericalIntegration::erf1c(double x) const {
  const static double alpha_erf1 = 1.0 / 0.30;
  return 0.5 * std::erfc(std::abs(x / (1.0 - x * x)) * alpha_erf1);
}

}  // namespace xtp
}  // namespace votca
