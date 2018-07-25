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



#include <votca/xtp/numerical_integrations.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/aoshell.h>
#include <votca/tools/constants.h>

#include <votca/xtp/aomatrix.h>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <iterator>
#include <string>




namespace votca {
    namespace xtp {
     

    double NumericalIntegration::getExactExchange(const std::string _functional) {      

      double exactexchange = 0.0;
      Vxc_Functionals map;
      std::vector<std::string> strs;

      boost::split(strs, _functional, boost::is_any_of(" "));
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
          fprintf(stderr, "Functional '%d' not found\n", func_id);
          exit(1);
        }
        if (exactexchange > 0 && func.cam_alpha > 0) {
          throw std::runtime_error("You have specified two functionals with exact exchange");
        }
        exactexchange += func.cam_alpha;
      }
      return exactexchange;

    }
    
    
   
        
  void NumericalIntegration::setXCfunctional(const std::string _functional) {

      Vxc_Functionals map;
      std::vector<std::string> strs;
      boost::split(strs, _functional, boost::is_any_of(" "));
      xfunc_id = 0;


      _use_separate = false;
      cfunc_id = 0;
      if (strs.size() == 1) {
        xfunc_id = map.getID(strs[0]);
      } else if (strs.size() == 2) {
        cfunc_id = map.getID(strs[0]);
        xfunc_id = map.getID(strs[1]);
        _use_separate = true;
      } else {
        std::cout << "LIBXC " << strs.size() << std::endl;
        throw std::runtime_error("LIBXC. Please specify one combined or an exchange and a correlation functionals");
      }
      
        if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
          fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
          exit(1);
        }
        xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
        if (xfunc.info->kind != 2 && !_use_separate) {
          throw std::runtime_error("Your functional misses either correlation or exchange, please specify another functional, separated by whitespace");
        }
        if (_use_separate) {
          if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
            fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
            exit(1);
          }
          xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED);
          xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
          if ((xfunc.info->kind + cfunc.info->kind) != 1) {
            throw std::runtime_error("Your functionals are not one exchange and one correlation");
          }
        }
      

      setXC = true;
      return;
    }

        
        void NumericalIntegration::EvaluateXC(const double rho, const Eigen::Vector3d& grad_rho, double& f_xc, double& df_drho, double& df_dsigma) {

        double sigma = (grad_rho.transpose()*grad_rho).value();
        double exc[1];
        double vsigma[1]; // libxc 
        double vrho[1]; // libxc df/drho
        switch (xfunc.info->family) {
          case XC_FAMILY_LDA:
            xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
            break;
          case XC_FAMILY_GGA:
          case XC_FAMILY_HYB_GGA:
            xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
            break;
        }
        f_xc = exc[0];
        df_drho = vrho[0];
        df_dsigma = vsigma[0];
        if (_use_separate) {
          // via libxc correlation part only
          switch (cfunc.info->family) {
            case XC_FAMILY_LDA:
              xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
              break;
            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
              xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
              break;
          }

          f_xc += exc[0];
          df_drho += vrho[0];
          df_dsigma += vsigma[0];
        }
  
      return;
    }
        
        
    double NumericalIntegration::IntegratePotential(const tools::vec& rvector) {

      double result = 0.0;
      assert(density_set && "Density not calculated");
      for (unsigned i = 0; i < _grid_boxes.size(); i++) {
        const std::vector<tools::vec>& points = _grid_boxes[i].getGridPoints();
        const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
        const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
        for (unsigned j = 0; j < points.size(); j++) {
          double dist = abs(points[j] - rvector);
          result -= weights[j] * densities[j] / dist;
        }
      }
      return result;
    }
               
        
        
        
  void NumericalIntegration::SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid) {
      const double boxsize = 1;//1 bohr

      std::vector< std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > > boxes;
      tools::vec min = tools::vec(std::numeric_limits<double>::max());
      tools::vec max = tools::vec(std::numeric_limits<double>::min());

      for (unsigned i = 0; i < grid.size(); i++) {
        for (unsigned j = 0; j < grid[i].size(); j++) {
          const tools::vec& pos = grid[i][j].grid_pos;
          if (pos.getX() > max.getX()) {
            max.x() = pos.getX();
          } else if (pos.getX() < min.getX()) {
            min.x() = pos.getX();
          }
          if (pos.getY() > max.getY()) {
            max.y() = pos.getY();
          } else if (pos.getY() < min.getY()) {
            min.y() = pos.getY();
          }
          if (pos.getZ() > max.getZ()) {
            max.z() = pos.getZ();
          } else if (pos.getZ() < min.getZ()) {
            min.z() = pos.getZ();
          }
        }
      }

      tools::vec molextension = (max - min);
      tools::vec numberofboxes = molextension / boxsize;
      tools::vec roundednumofbox = tools::vec(std::ceil(numberofboxes.getX()), std::ceil(numberofboxes.getY()), std::ceil(numberofboxes.getZ()));

      //creating temparray
      for (unsigned i = 0; i<unsigned(roundednumofbox.getX()); i++) {
        std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > boxes_yz;
        for (unsigned j = 0; j<unsigned(roundednumofbox.getY()); j++) {
          std::vector< std::vector< GridContainers::integration_grid* > > boxes_z;
          for (unsigned k = 0; k<unsigned(roundednumofbox.getZ()); k++) {
            std::vector< GridContainers::integration_grid* > box;
            box.reserve(100);
            boxes_z.push_back(box);
          }
          boxes_yz.push_back(boxes_z);
        }
        boxes.push_back(boxes_yz);
      }

      for (auto & atomgrid : grid) {
        for (auto & gridpoint : atomgrid) {
          tools::vec pos = gridpoint.grid_pos - min;
          tools::vec index = pos / boxsize;
          int i_x = int(index.getX());
          int i_y = int(index.getY());
          int i_z = int(index.getZ());
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

            for (const auto&point : box) {
              gridbox.addGridPoint(*point);
            }
            _grid_boxes.push_back(gridbox);
          }
        }
      }
      return;
    }
        
  void NumericalIntegration::FindSignificantShells() {

      for (unsigned i = 0; i < _grid_boxes.size(); ++i) {
        GridBox & box = _grid_boxes[i];
        for (const AOShell* store:(*_basis)) {
          const double decay = store->getMinDecay();
          const tools::vec& shellpos =store->getPos();
          for (const auto& point : box.getGridPoints()) {
            tools::vec dist = shellpos - point;
            double distsq = dist*dist;
            // if contribution is smaller than -ln(1e-10), add atom to list
            if ((decay * distsq) < 20.7) {
              box.addShell(store);
              break;
            }
          }
        }
      }

      std::vector< GridBox > _grid_boxes_copy;

      int combined = 0;
      std::vector<bool> Compared = std::vector<bool>(_grid_boxes.size(), false);
      for (unsigned i = 0; i < _grid_boxes.size(); i++) {
        if (Compared[i]) {
          continue;
        }
        GridBox box = _grid_boxes[i];
        if (box.Shellsize() < 1) {
          continue;
        }
        Compared[i] = true;
        for (unsigned j = i + 1; j < _grid_boxes.size(); j++) {
          if (GridBox::compareGridboxes(_grid_boxes[i], _grid_boxes[j])) {
            Compared[j] = true;
            box.addGridBox(_grid_boxes[j]);
            combined++;
          }

        }
        _grid_boxes_copy.push_back(box);
      }

      std::vector<unsigned> sizes;
      sizes.reserve(_grid_boxes_copy.size());
      for (auto& box : _grid_boxes_copy) {
        sizes.push_back(box.size() * box.Matrixsize());
      }
      std::vector<unsigned> indexes = std::vector<unsigned>(sizes.size());
      std::iota(indexes.begin(), indexes.end(), 0);
      std::sort(indexes.begin(), indexes.end(), [&sizes](unsigned i1, unsigned i2) {
        return sizes[i1] > sizes[i2];
      });

      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<unsigned> scores = std::vector<unsigned>(nthreads, 0);
      std::vector< std::vector<unsigned> > indices;
      for (unsigned i = 0; i < nthreads; ++i) {
        std::vector<unsigned> thread_box_indices;
        indices.push_back(thread_box_indices);
      }
      for (const auto index : indexes) {
        unsigned thread = 0;
        unsigned minimum = std::numeric_limits<unsigned>::max();
        for (unsigned i = 0; i < scores.size(); ++i) {
          if (scores[i] < minimum) {
            minimum = scores[i];
            thread = i;
          }
        }
        indices[thread].push_back(index);
        scores[thread] += sizes[index];
      }
      thread_start = std::vector<unsigned>(0);
      thread_stop = std::vector<unsigned>(0);
      unsigned start = 0;
      unsigned stop = 0;
      unsigned indexoffirstgridpoint = 0;
      _grid_boxes.resize(0);
      for (const std::vector<unsigned>& thread_index : indices) {
        thread_start.push_back(start);
        stop = start + thread_index.size();
        thread_stop.push_back(stop);
        start = stop;
        for (const unsigned index : thread_index) {
          GridBox newbox = _grid_boxes_copy[index];
          newbox.setIndexoffirstgridpoint(indexoffirstgridpoint);
          indexoffirstgridpoint += newbox.size();
          newbox.PrepareForIntegration();
          _grid_boxes.push_back(newbox);
        }
      }
      return;
    }
        
        
        
  Eigen::MatrixXd NumericalIntegration::IntegrateVXC(const Eigen::MatrixXd& _density_matrix) {
      Eigen::MatrixXd Vxc = Eigen::MatrixXd::Zero(_density_matrix.rows(), _density_matrix.cols());
      EXC = 0;
      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >vxc_thread;
      std::vector<double> Exc_thread = std::vector<double>(nthreads, 0.0);
      for (unsigned i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd Vxc_thread = Eigen::MatrixXd::Zero(_density_matrix.rows(), _density_matrix.cols());
        vxc_thread.push_back(Vxc_thread);
      }
      

#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {

          double EXC_box = 0.0;
          GridBox& box = _grid_boxes[i];
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(_density_matrix);
          const Eigen::MatrixXd DMAT_symm= DMAT_here+DMAT_here.transpose();
          double cutoff=1.e-40/_density_matrix.rows()/_density_matrix.rows();
          if (DMAT_here.cwiseAbs2().maxCoeff()<cutoff ){
            continue;
          }
          
          Eigen::MatrixXd Vxc_here = Eigen::MatrixXd::Zero(DMAT_here.rows(), DMAT_here.cols());
          const std::vector<tools::vec>& points = box.getGridPoints();
          const std::vector<double>& weights = box.getGridWeights();

          //iterate over gridpoints
          for (unsigned p = 0; p < box.size(); p++) {
            Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
            Eigen::MatrixX3d ao_grad= Eigen::MatrixX3d::Zero(box.Matrixsize(),3);
            const std::vector<GridboxRange>& aoranges = box.getAOranges();
            const std::vector<const AOShell* >& shells = box.getShells();
                        
            for (unsigned j = 0; j < box.Shellsize(); ++j) {
              Eigen::Block<Eigen::MatrixX3d> grad_block=ao_grad.block(aoranges[j].start,0,aoranges[j].size,3);
              Eigen::VectorBlock<Eigen::VectorXd> ao_block=ao.segment(aoranges[j].start,aoranges[j].size);
              shells[j]->EvalAOspace(ao_block,grad_block,points[p]);             
            }
            double rho = (ao.transpose()*DMAT_here*ao).value();
            Eigen::Vector3d rho_grad = ao.transpose()*DMAT_symm*ao_grad;
            double weight = weights[p];
            if (rho*weight < 1.e-20) continue; // skip the rest, if density is very small
            double f_xc; // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
            double df_drho; // v_xc_rho(r) = df/drho
            double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
            EvaluateXC(rho, rho_grad, f_xc, df_drho, df_dsigma);          
            auto _addXC = weight * (0.5*df_drho * ao+ 2.0 * df_dsigma *ao_grad*rho_grad);
            // Exchange correlation energy
            EXC_box += weight * rho * f_xc;
            Vxc_here.noalias() += _addXC* ao.transpose();
          }
          box.AddtoBigMatrix(vxc_thread[thread], Vxc_here);
          Exc_thread[thread] += EXC_box;
        }
      }
      for (unsigned i = 0; i < nthreads; ++i) {
        Vxc += vxc_thread[i];
        EXC += Exc_thread[i];
      }
    
      return Vxc+Vxc.transpose();
    }
  
  
   Eigen::MatrixXd NumericalIntegration::IntegrateExternalPotential(const std::vector<double>& Potentialvalues) {

      Eigen::MatrixXd ExternalMat = Eigen::MatrixXd::Zero(_basis->AOBasisSize(), _basis->AOBasisSize());
      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >vex_thread;
      std::vector<double> Exc_thread = std::vector<double>(nthreads, 0.0);
      for (unsigned i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd Vex_thread = Eigen::MatrixXd::Zero(ExternalMat.rows(), ExternalMat.cols());
        vex_thread.push_back(Vex_thread);
      }


#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {

          GridBox& box = _grid_boxes[i];
          Eigen::MatrixXd Vex_here = Eigen::MatrixXd::Zero(box.Matrixsize(), box.Matrixsize());
          const std::vector<tools::vec>& points = box.getGridPoints();
          const std::vector<double>& weights = box.getGridWeights();

          //iterate over gridpoints
          for (unsigned p = 0; p < box.size(); p++) {
            Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
            const std::vector<GridboxRange>& aoranges = box.getAOranges();
            const std::vector<const AOShell* > shells = box.getShells();
            for (unsigned j = 0; j < box.Shellsize(); ++j) {
              Eigen::VectorBlock<Eigen::VectorXd> ao_block=ao.segment(aoranges[j].start,aoranges[j].size);
              shells[j]->EvalAOspace(ao_block, points[p]);
            }
            Eigen::VectorXd _addEX = weights[p] * Potentialvalues[box.getIndexoffirstgridpoint() + p] * ao;
            Vex_here += _addEX.transpose() * ao;
          }
          box.AddtoBigMatrix(vex_thread[thread], Vex_here);

        }
      }
      for (unsigned i = 0; i < nthreads; ++i) {
        ExternalMat += vex_thread[i];
      }

      ExternalMat += ExternalMat.transpose();
      return ExternalMat;

    }
      
  double NumericalIntegration::IntegrateDensity(const Eigen::MatrixXd& _density_matrix) {
      double N = 0;
      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);

#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {

          double N_box = 0.0;
          GridBox& box = _grid_boxes[i];
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(_density_matrix);
          const std::vector<tools::vec>& points = box.getGridPoints();
          const std::vector<double>& weights = box.getGridWeights();
          box.prepareDensity();

          //iterate over gridpoints
          for (unsigned p = 0; p < box.size(); p++) {
            Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
            const std::vector<GridboxRange>& aoranges = box.getAOranges();
            const std::vector<const AOShell* > shells = box.getShells();
            for (unsigned j = 0; j < box.Shellsize(); ++j) {
              Eigen::VectorBlock<Eigen::VectorXd> ao_block=ao.segment(aoranges[j].start,aoranges[j].size);
              shells[j]->EvalAOspace(ao_block, points[p]);
            }
            double rho =(ao.transpose()*DMAT_here*ao)(0, 0);
            box.addDensity(rho);
            N_box += rho * weights[p];
          }
          N_thread[thread] += N_box;
        }
      }
      for (unsigned i = 0; i < nthreads; ++i) {
        N += N_thread[i];
      }
      density_set = true;
      return N;
    }
        
      Gyrationtensor NumericalIntegration::IntegrateGyrationTensor(const Eigen::MatrixXd& _density_matrix) {

      double N = 0;
      tools::vec centroid = tools::vec(0.0);
      tools::matrix gyration = tools::matrix(0.0);
      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);
      // centroid
      std::vector<tools::vec> centroid_thread;
      std::vector<tools::matrix> gyration_thread;
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        tools::vec tempvec = tools::vec(0.0);
        centroid_thread.push_back(tempvec);
        tools::matrix tempmatrix = tools::matrix(0.0);
        gyration_thread.push_back(tempmatrix);
      }

#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {
          double N_box = 0.0;
          tools::vec centroid_box = tools::vec(0.0);
          tools::matrix gyration_box = tools::matrix(0.0);
          GridBox& box = _grid_boxes[i];
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(_density_matrix);
          const std::vector<tools::vec>& points = box.getGridPoints();
          const std::vector<double>& weights = box.getGridWeights();
          box.prepareDensity();
          //iterate over gridpoints
          for (unsigned p = 0; p < box.size(); p++) {
            Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
            const std::vector<GridboxRange>& aoranges = box.getAOranges();
            const std::vector<const AOShell* > shells = box.getShells();
            for (unsigned j = 0; j < box.Shellsize(); ++j) {
              Eigen::VectorBlock<Eigen::VectorXd> ao_block=ao.segment(aoranges[j].start,aoranges[j].size);
              shells[j]->EvalAOspace(ao_block, points[p]);
            }
            double rho =(ao.transpose()*DMAT_here*ao)(0, 0);
            box.addDensity(rho);
            N_box += rho * weights[p];
            centroid_box+=rho * weights[p] * points[p];
            gyration_box += rho * weights[p] * (points[p]|points[p]); 
          }
          N_thread[thread] += N_box;
          centroid_thread[thread] += centroid_box;
          gyration_thread[thread] += gyration_box;
        }
      }
      for (unsigned i = 0; i < nthreads; ++i) {
        N += N_thread[i];
        centroid += centroid_thread[i];
        gyration += gyration_thread[i];
      }
      density_set = true;

      // Normalize
      centroid = centroid / N;
      gyration = gyration / N;
      gyration=gyration-(centroid|centroid);
      Gyrationtensor gyro;
      gyro.mass=N;
      gyro.centroid=centroid;
      gyro.gyration=gyration;
      
      return gyro;
    }

        
std::vector<const tools::vec *> NumericalIntegration::getGridpoints() const{
    std::vector<const tools::vec *> gridpoints;
    for (unsigned i = 0; i < _grid_boxes.size(); i++) {
      const std::vector<tools::vec>& points = _grid_boxes[i].getGridPoints();
      for (unsigned j = 0; j < points.size(); j++) {
        gridpoints.push_back(&points[j]);

      }
    }
    return gridpoints;
  }


Eigen::MatrixXd NumericalIntegration::IntegratePotential(const AOBasis& externalbasis){
  Eigen::MatrixXd Potential=Eigen::MatrixXd::Zero(externalbasis.AOBasisSize(),externalbasis.AOBasisSize());
  
  assert(density_set && "Density not calculated");
  for (unsigned i = 0; i < _grid_boxes.size(); i++) {
    const std::vector<tools::vec>& points = _grid_boxes[i].getGridPoints();
    const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
    const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
    for (unsigned j = 0; j < points.size(); j++) {
      double weighteddensity=weights[j]*densities[j];
      if (weighteddensity<1e-12){
        continue;
      }
      AOESP esp;
      esp.setPosition(points[j]);
      esp.Fill(externalbasis);
      Potential+=weighteddensity*esp.Matrix();
    }
  }
  return Potential; 
}
        
        
void NumericalIntegration::GridSetup(std::string type, std::vector<QMAtom*> _atoms,const AOBasis* basis) {
      _basis = basis;
      std::vector< std::vector< GridContainers::integration_grid > > grid;
      const double pi = boost::math::constants::pi<double>();
      // get GridContainer
      GridContainers initialgrids;
      // get radial grid per element
      EulerMaclaurinGrid _radialgrid;
      _radialgrid.getRadialGrid(basis, _atoms, type, initialgrids); // this checks out 1:1 with NWChem results! AWESOME
      std::map<std::string, GridContainers::radial_grid>::iterator it;
      LebedevGrid _sphericalgrid;
      for (it = initialgrids._radial_grids.begin(); it != initialgrids._radial_grids.end(); ++it) {
        _sphericalgrid.getSphericalGrid(_atoms, type, initialgrids);
      }
      // for the partitioning, we need all inter-center distances later, stored in one-directional list
      int ij = 0;
      Rij.push_back(0.0); // 1st center "self-distance"
      std::vector< QMAtom* > ::iterator ait;
      std::vector< QMAtom* > ::iterator bit;
      int i = 1;
      for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
        // get center coordinates in Bohr
        tools::vec pos_a = (*ait)->getPos();
        int j = 0;
        for (bit = _atoms.begin(); bit != ait; ++bit) {
          ij++;
          // get center coordinates in Bohr
          tools::vec pos_b = (*bit)->getPos();
          Rij.push_back(1.0 / abs(pos_a - pos_b));
          j++;
        } // atoms
        Rij.push_back(0.0); // self-distance again
        i++;
      } // atoms

      int i_atom = 0;
      _totalgridsize = 0;
      for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
        // get center coordinates in Bohr
        std::vector< GridContainers::integration_grid > _atomgrid;
        const tools::vec & atomA_pos = (*ait)->getPos();
        const std::string & name = (*ait)->getType();
        // get radial grid information for this atom type
        GridContainers::radial_grid _radial_grid = initialgrids._radial_grids.at(name);
        // get spherical grid information for this atom type
        GridContainers::spherical_grid _spherical_grid = initialgrids._spherical_grids.at(name);
        // maximum order (= number of points) in spherical integration grid
        int maxorder = _sphericalgrid.Type2MaxOrder(name, type);
        int maxindex = _sphericalgrid.getIndexFromOrder(maxorder);
        // for pruning of integration grid, get interval boundaries for this element
        std::vector<double> PruningIntervals = _radialgrid.getPruningIntervals(name);
        int current_order = 0;
        // get spherical grid
        std::vector<double> _theta;
        std::vector<double> _phi;
        std::vector<double> _weight;
        // for each radial value
        for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
          double r = _radial_grid.radius[_i_rad];
          int order;
          // which Lebedev order for this point?
          if (maxindex == 1) {
            // smallest possible grid anyway, nothing to do
            order = maxorder;
          } else if (maxindex == 2) {
            // only three intervals
            if (r < PruningIntervals[0]) {
              order = _sphericalgrid.getOrderFromIndex(1); //1;
            } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[3])) {
              order = _sphericalgrid.getOrderFromIndex(2);
            } else {
              order = _sphericalgrid.getOrderFromIndex(1);
            } // maxorder == 2
          } else {
            // five intervals
            if (r < PruningIntervals[0]) {
              order = _sphericalgrid.getOrderFromIndex(int(2));
            } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[1])) {
              order = _sphericalgrid.getOrderFromIndex(4);
            } else if ((r >= PruningIntervals[1]) && (r < PruningIntervals[2])) {
              order = _sphericalgrid.getOrderFromIndex(std::max(maxindex - 1, 4));
            } else if ((r >= PruningIntervals[2]) && (r < PruningIntervals[3])) {
              order = maxorder;
            } else {
              order = _sphericalgrid.getOrderFromIndex(std::max(maxindex - 1, 1));
            }
          }

          // get new spherical grid, if order changed
          if (order != current_order) {
            _theta.clear();
            _phi.clear();
            _weight.clear();
            _sphericalgrid.getUnitSphereGrid(order, _theta, _phi, _weight);
            current_order = order;
          }

          for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++) {
            double p = _phi[_i_sph] * pi / 180.0; // back to rad
            double t = _theta[_i_sph] * pi / 180.0; // back to rad
            double ws = _weight[_i_sph];
            const tools::vec s = tools::vec(sin(p) * cos(t), sin(p) * sin(t), cos(p));
            GridContainers::integration_grid _gridpoint;
            _gridpoint.grid_pos = atomA_pos + r*s;
            _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;
            _atomgrid.push_back(_gridpoint);
          } // spherical gridpoints
        } // radial gridpoint


        // get all distances from grid points to centers
        std::vector< std::vector<double> > rq;
        // for each center
        for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {
          // get center coordinates
          const tools::vec & atom_pos = (*bit)->getPos();
          std::vector<double> temp;
          // for each gridpoint
          for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {
            temp.push_back(abs(git->grid_pos - atom_pos));
          } // gridpoint of _atomgrid
          rq.push_back(temp);
        } // centers
        
        // find nearest-neighbor of this atom
        double distNN = std::numeric_limits<double>::max();
        std::vector< QMAtom* > ::iterator NNit;
              // now check all other centers
        int i_b = 0;
        for (bit = _atoms.begin(); bit != _atoms.end(); ++bit) {
          if (bit != ait) {
            // get center coordinates
            const tools::vec & atomB_pos = (*bit)->getPos();
            double distSQ = (atomA_pos - atomB_pos)*(atomA_pos - atomB_pos);
            // update NN distance and iterator
            if (distSQ < distNN) {
              distNN = distSQ;
              NNit = bit;
            }
          } // if ( ait != bit) 
          i_b++;
        }// bit centers
#pragma omp parallel for schedule(guided)
        for (unsigned i_grid = 0; i_grid < _atomgrid.size(); i_grid++) {
          // call some shit called grid_ssw0 in NWChem
          std::vector<double> _p = SSWpartition(i_grid, _atoms.size(), rq);
          // check weight sum
          double wsum = 0.0;
          for (const auto&p:_p) {
            wsum +=p;
          }
          if (wsum != 0.0) {
            // update the weight of this grid point
            _atomgrid[i_grid].grid_weight *= _p[i_atom] / wsum;
          } else {
            std::cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
            throw std::runtime_error("\nThis should never happen!");
          }
        } // partition weight for each gridpoint

        // now remove points from the grid with negligible weights
        for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();) {
          if (git->grid_weight < 1e-13) {
            git = _atomgrid.erase(git);
          } else {
            ++git;
          }
        }
        _totalgridsize += _atomgrid.size();
        grid.push_back(_atomgrid);
        i_atom++;
      } // atoms
      SortGridpointsintoBlocks(grid);
      FindSignificantShells();
      return;
    }

    std::vector<double> NumericalIntegration::SSWpartition(int igrid, int ncenters,const std::vector< std::vector<double> >& rq) {
      const double ass = 0.725;
      // initialize partition vector to 1.0
      std::vector<double> p(ncenters, 1.0);
      const double tol_scr = 1e-10;
      const double leps = 1e-6;
      // go through centers
      for (int i = 1; i < ncenters; i++) {
        int ij = i * (i + 1) / 2 - 1; // indexing magic
        double rag = rq[i][igrid];
        // through all other centers (one-directional)
        for (int j = 0; j < i; j++) {
          ij++;
          if ((std::abs(p[i]) > tol_scr) || (std::abs(p[j]) > tol_scr)) {
            double mu = (rag - rq[j][igrid]) * Rij[ij];
            if (mu > ass) {
              p[i] = 0.0;
            } else if (mu < -ass) {
              p[j] = 0.0;
            } else {
              double sk;
              if (std::abs(mu) < leps) {
                sk = -1.88603178008 * mu + 0.5;
              } else {
                sk =erf1c(mu);
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

    double NumericalIntegration::erf1c(double x){ 
        const static double alpha_erf1=1.0/0.30;
        return 0.5*std::erfc(std::abs(x/(1.0-x*x))*alpha_erf1);              
    }
              
   
                                                                                                
    }
}
