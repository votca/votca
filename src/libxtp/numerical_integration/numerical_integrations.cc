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

    NumericalIntegration::~NumericalIntegration(){
      if (_setXC) {
        xc_func_end(&xfunc);
        if (_use_separate) {
          xc_func_end(&cfunc);
        }
      }
    }

    

    double NumericalIntegration::getExactExchange(const std::string& functional) {      

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
         throw std::runtime_error((boost::format("Functional %s not found\n") %strs[i]).str());
        }
        if (exactexchange > 0 && func.cam_alpha > 0) {
          throw std::runtime_error("You have specified two functionals with exact exchange");
        }
        exactexchange += func.cam_alpha;
        xc_func_end(&func);
      }
      
      return exactexchange;

    }
    
    
   void NumericalIntegration::setXCfunctional(const std::string& functional) {

      Vxc_Functionals map;
      std::vector<std::string> strs;
      tools::Tokenizer tok(functional," ,\n\t");
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
        throw std::runtime_error("LIBXC. Please specify one combined or an exchange and a correlation functionals");
      }

      if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
        throw std::runtime_error((boost::format("Functional %s not found\n") %strs[0]).str());
      }
      if (xfunc.info->kind != 2 && !_use_separate) {
        throw std::runtime_error("Your functional misses either correlation or exchange, please specify another functional, separated by whitespace");
      }
      if (_use_separate) {
        if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
          throw std::runtime_error((boost::format("Functional %s not found\n") %strs[1]).str());
        }
        if ((xfunc.info->kind + cfunc.info->kind) != 1) {
          throw std::runtime_error("Your functionals are not one exchange and one correlation");
        }
      }
      _setXC = true;
      return;
    }

        
        void NumericalIntegration::EvaluateXC(const double rho, const double sigma, double& f_xc, double& df_drho, double& df_dsigma) {

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
      assert(_density_set && "Density not calculated");
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
               
        
        
        
  void NumericalIntegration::SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::Cartesian_gridpoint > >& grid) {
      const double boxsize = 1;//1 bohr

      std::vector< std::vector< std::vector< std::vector< GridContainers::Cartesian_gridpoint* > > > > boxes;
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
        std::vector< std::vector< std::vector< GridContainers::Cartesian_gridpoint* > > > boxes_yz;
        for (unsigned j = 0; j<unsigned(roundednumofbox.getY()); j++) {
          std::vector< std::vector< GridContainers::Cartesian_gridpoint* > > boxes_z;
          for (unsigned k = 0; k<unsigned(roundednumofbox.getZ()); k++) {
            std::vector< GridContainers::Cartesian_gridpoint* > box;
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
        
  void NumericalIntegration::FindSignificantShells(const AOBasis& basis) {
      for (unsigned i = 0; i < _grid_boxes.size(); ++i) {
        GridBox & box = _grid_boxes[i];
        for (const AOShell* store:basis) {
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

      std::vector< GridBox > grid_boxes_copy;
      int combined = 0;
      //use vecot of bool to indicate if a gridbox has already been merged into another
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
      std::vector<unsigned> sizes;
      sizes.reserve(grid_boxes_copy.size());
      for (auto& box : grid_boxes_copy) {
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
          GridBox newbox = grid_boxes_copy[index];
          newbox.setIndexoffirstgridpoint(indexoffirstgridpoint);
          indexoffirstgridpoint += newbox.size();
          newbox.PrepareForIntegration();
          _grid_boxes.push_back(newbox);
        }
      }
      return;
    }
        
        
        
  Eigen::MatrixXd NumericalIntegration::IntegrateVXC(const Eigen::MatrixXd& density_matrix) {
      Eigen::MatrixXd Vxc = Eigen::MatrixXd::Zero(density_matrix.rows(), density_matrix.cols());
      _EXC = 0;
      unsigned nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      std::vector<Eigen::MatrixXd >vxc_thread;
      std::vector<double> Exc_thread = std::vector<double>(nthreads, 0.0);
      for (unsigned i = 0; i < nthreads; ++i) {
        Eigen::MatrixXd Vxc_thread = Eigen::MatrixXd::Zero(density_matrix.rows(), density_matrix.cols());
        vxc_thread.push_back(Vxc_thread);
      }
      
#pragma omp parallel for
      for (unsigned thread = 0; thread < nthreads; ++thread) {
        for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {

          double EXC_box = 0.0;
          const GridBox& box = _grid_boxes[i];
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
          const Eigen::MatrixXd DMAT_symm= DMAT_here+DMAT_here.transpose();
          double cutoff=1.e-40/density_matrix.rows()/density_matrix.rows();
          if (DMAT_here.cwiseAbs2().maxCoeff()<cutoff ){
            continue;
          }     
          Eigen::MatrixXd Vxc_here = Eigen::MatrixXd::Zero(DMAT_here.rows(), DMAT_here.cols());
          const std::vector<tools::vec>& points = box.getGridPoints();
          const std::vector<double>& weights = box.getGridWeights();
          Eigen::VectorXd ao = Eigen::VectorXd::Zero(box.Matrixsize());
          Eigen::MatrixX3d ao_grad= Eigen::MatrixX3d::Zero(box.Matrixsize(),3);
          const std::vector<GridboxRange>& aoranges = box.getAOranges();
          const std::vector<const AOShell* >& shells = box.getShells();
          Eigen::VectorXd grad =Eigen::VectorXd::Zero(box.Matrixsize());
          //iterate over gridpoints
          for (unsigned p = 0; p < box.size(); p++) {
            ao.setZero(box.Matrixsize());
            ao_grad.setZero(box.Matrixsize(),3);
            for (unsigned j = 0; j < box.Shellsize(); ++j) {
              Eigen::Block<Eigen::MatrixX3d> grad_block=ao_grad.block(aoranges[j].start,0,aoranges[j].size,3);
              Eigen::VectorBlock<Eigen::VectorXd> ao_block=ao.segment(aoranges[j].start,aoranges[j].size);
              shells[j]->EvalAOspace(ao_block,grad_block,points[p]);             
            }
            const double rho =0.5* (ao.transpose()*DMAT_symm*ao).value();
            const double weight = weights[p];
            if (rho*weight < 1.e-20) continue; // skip the rest, if density is very small
            const Eigen::Vector3d rho_grad = ao.transpose()*DMAT_symm*ao_grad;
            grad=ao_grad*rho_grad;
            const double sigma = (rho_grad.transpose()*rho_grad).value();
            double f_xc; // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
            double df_drho; // v_xc_rho(r) = df/drho
            double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
            EvaluateXC(rho, sigma, f_xc, df_drho, df_dsigma); 
            EXC_box += weight * rho * f_xc;
            auto addXC = weight * (0.5*df_drho * ao+ 2.0 * df_dsigma *grad);
            // Exchange correlation energy
            Vxc_here.noalias() += addXC* ao.transpose();
          }
          box.AddtoBigMatrix(vxc_thread[thread], Vxc_here);
          Exc_thread[thread] += EXC_box;
        }
      }
      for (unsigned i = 0; i < nthreads; ++i) {
        Vxc += vxc_thread[i];
        _EXC += Exc_thread[i];
      }    
      return Vxc+Vxc.transpose();
    }
  
  
   Eigen::MatrixXd NumericalIntegration::IntegrateExternalPotential(const std::vector<double>& Potentialvalues) {

      Eigen::MatrixXd ExternalMat = Eigen::MatrixXd::Zero(_AOBasisSize, _AOBasisSize);
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

          const GridBox& box = _grid_boxes[i];
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
            Eigen::VectorXd addEX = weights[p] * Potentialvalues[box.getIndexoffirstgridpoint() + p] * ao;
            Vex_here += addEX.transpose() * ao;
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
      
  double NumericalIntegration::IntegrateDensity(const Eigen::MatrixXd& density_matrix) {
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
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
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
      _density_set = true;
      return N;
    }
        
      Gyrationtensor NumericalIntegration::IntegrateGyrationTensor(const Eigen::MatrixXd& density_matrix) {
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
          const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
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
      _density_set = true;
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
  
  assert(_density_set && "Density not calculated");
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
        

Eigen::MatrixXd NumericalIntegration::CalcInverseAtomDist(std::vector<QMAtom*>& atoms){
  Eigen::MatrixXd result=Eigen::MatrixXd::Zero(atoms.size(),atoms.size());
#pragma omp parallel for
  for (unsigned i=0;i<atoms.size();++i) {
    QMAtom* atom_a=atoms[i];
    const tools::vec& pos_a = atom_a->getPos();
   for (unsigned j=0;j<i;++j) {
      QMAtom* atom_b=atoms[j];
      const tools::vec& pos_b = atom_b->getPos();
      result(j,i)=1/tools::abs(pos_a-pos_b);
    } 
  }
  return result+result.transpose();
}

int NumericalIntegration::UpdateOrder(LebedevGrid& sphericalgridofElement, int maxorder, std::vector<double>& PruningIntervals, double r){
  int order;
  int maxindex=sphericalgridofElement.getIndexFromOrder(maxorder);
  if (maxindex == 1) {
    // smallest possible grid anyway, nothing to do
    order = maxorder;
  } else if (maxindex == 2) {
    // only three intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(1); //1;
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[3])) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else {
      order = sphericalgridofElement.getOrderFromIndex(1);
    } // maxorder == 2
  } else {
    // five intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[1])) {
      order = sphericalgridofElement.getOrderFromIndex(4);
    } else if ((r >= PruningIntervals[1]) && (r < PruningIntervals[2])) {
      order = sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 4));
    } else if ((r >= PruningIntervals[2]) && (r < PruningIntervals[3])) {
      order = maxorder;
    } else {
      order = sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 1));
    }
  }
  return order;
    }

    GridContainers::Cartesian_gridpoint NumericalIntegration::CreateCartesianGridpoint(const tools::vec& atomA_pos,
            GridContainers::radial_grid& radial_grid, GridContainers::spherical_grid& spherical_grid,
            unsigned i_rad, unsigned i_sph) {
      GridContainers::Cartesian_gridpoint gridpoint;
      double p = spherical_grid.phi[i_sph];
      double t = spherical_grid.theta[i_sph];
      const tools::vec s = tools::vec(sin(p) * cos(t), sin(p) * sin(t), cos(p));
      double r = radial_grid.radius[i_rad];
      gridpoint.grid_pos = atomA_pos + r*s;
      gridpoint.grid_weight = radial_grid.weight[i_rad] * spherical_grid.weight[i_sph];
      return gridpoint;
    }

    Eigen::MatrixXd NumericalIntegration::CalcDistanceAtomsGridpoints(std::vector<QMAtom*>& atoms, std::vector<GridContainers::Cartesian_gridpoint>& atomgrid){
      Eigen::MatrixXd result=Eigen::MatrixXd::Zero(atoms.size(),atomgrid.size());
     #pragma omp parallel for
      for (unsigned i=0;i<atoms.size();++i) {
        QMAtom* atom=atoms[i];
        const tools::vec & atom_pos =atom->getPos();
        for (unsigned j=0;j<atomgrid.size();++j) {
          const auto& gridpoint=atomgrid[j];
          result(i,j)=tools::abs(atom_pos-gridpoint.grid_pos);
        } 
      }
      return result;
    }

    void NumericalIntegration::SSWpartitionAtom(std::vector<QMAtom*>& atoms, std::vector<GridContainers::Cartesian_gridpoint>& atomgrid
                                                , unsigned i_atom, const Eigen::MatrixXd& Rij){
      Eigen::MatrixXd AtomGridDist=CalcDistanceAtomsGridpoints(atoms, atomgrid);
      
#pragma omp parallel for schedule(guided)
      for (unsigned i_grid = 0; i_grid < atomgrid.size(); i_grid++) {
        Eigen::VectorXd p = SSWpartition(i_grid, AtomGridDist,Rij);
        // check weight sum
        double wsum = p.sum();
        if (wsum != 0.0) {
          // update the weight of this grid point
          atomgrid[i_grid].grid_weight *= p[i_atom] / wsum;
        } else {
          std::cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
          throw std::runtime_error("\nThis should never happen!");
        }
      } // partition weight for each gridpoint
    }
        
void NumericalIntegration::GridSetup(const std::string& type, std::vector<QMAtom*> atoms,const AOBasis& basis) {
      _AOBasisSize=basis.AOBasisSize();
      GridContainers initialgrids;
      // get radial grid per element
      EulerMaclaurinGrid radialgridofElement;
      initialgrids.radial_grids=radialgridofElement.CalculateAtomicRadialGrids(basis, atoms, type); // this checks out 1:1 with NWChem results! AWESOME
      LebedevGrid sphericalgridofElement;
      initialgrids.spherical_grids=sphericalgridofElement.CalculateSphericalGrids(atoms,type);
      
      // for the partitioning, we need all inter-center distances later, stored in matrix
      Eigen::MatrixXd Rij=CalcInverseAtomDist(atoms);
      _totalgridsize = 0;
      std::vector< std::vector< GridContainers::Cartesian_gridpoint > > grid;
      
      for (unsigned i_atom=0;i_atom<atoms.size();++i_atom) {
        QMAtom* atom=atoms[i_atom];

        const tools::vec & atomA_pos =atom->getPos();
        const std::string & name = atom->getType();
        GridContainers::radial_grid radial_grid = initialgrids.radial_grids.at(name);
        GridContainers::spherical_grid spherical_grid = initialgrids.spherical_grids.at(name);
                
        // maximum order (= number of points) in spherical integration grid
        int maxorder = sphericalgridofElement.Type2MaxOrder(name, type);
        // for pruning of integration grid, get interval boundaries for this element
        std::vector<double> PruningIntervals = radialgridofElement.CalculatePruningIntervals(name);
        int current_order = 0;
        // for each radial value
        std::vector< GridContainers::Cartesian_gridpoint > atomgrid;
        for (unsigned i_rad = 0; i_rad < radial_grid.radius.size(); i_rad++) {
          double r = radial_grid.radius[i_rad];

          // which Lebedev order for this point?
          int order=UpdateOrder(sphericalgridofElement, maxorder, PruningIntervals, r);
          // get new spherical grid, if order changed
          if (order != current_order) {
            spherical_grid=sphericalgridofElement.CalculateUnitSphereGrid(order);
            current_order = order;
          }

          for (unsigned i_sph = 0; i_sph < spherical_grid.phi.size(); i_sph++) {
            GridContainers::Cartesian_gridpoint gridpoint=CreateCartesianGridpoint(atomA_pos, radial_grid, spherical_grid, i_rad,i_sph);
            atomgrid.push_back(gridpoint);
          } // spherical gridpoints
        } // radial gridpoint

        SSWpartitionAtom(atoms, atomgrid, i_atom,Rij);

        // now remove points from the grid with negligible weights
        for (std::vector<GridContainers::Cartesian_gridpoint >::iterator git = atomgrid.begin(); git != atomgrid.end();) {
          if (git->grid_weight < 1e-13) {
            git = atomgrid.erase(git);
          } else {
            ++git;
          }
        }
        _totalgridsize += atomgrid.size();
        grid.push_back(atomgrid);
      } // atoms
      SortGridpointsintoBlocks(grid);
      FindSignificantShells(basis);
      return;
    }

    Eigen::VectorXd NumericalIntegration::SSWpartition(int igrid, const Eigen::MatrixXd & rq,const Eigen::MatrixXd& Rij) {
      const double ass = 0.725;
      // initialize partition vector to 1.0
      Eigen::VectorXd p=Eigen::VectorXd::Ones(rq.rows());
      const double tol_scr = 1e-10;
      const double leps = 1e-6;
      // go through centers
      for (int i = 1; i < rq.rows(); i++) {
        double rag = rq(i,igrid);
        // through all other centers (one-directional)
        for (int j = 0; j < i; j++) {
          if ((std::abs(p[i]) > tol_scr) || (std::abs(p[j]) > tol_scr)) {
            double mu = (rag - rq(j,igrid)) * Rij(j,i);
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
