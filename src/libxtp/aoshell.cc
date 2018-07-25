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
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include <votca/xtp/aomatrix.h>

namespace votca { namespace xtp {
    
  void AOShell::normalizeContraction(){
    
    AOOverlap overlap;
    Eigen::MatrixXd block=overlap.FillShell(this);
    std::vector<int> numsubshell=NumFuncSubShell(_type);
    int contraction_index=_Lmin;
    int aoindex=0;
    for (unsigned i=0;i<numsubshell.size();++i){
      double norm=std::sqrt(block(aoindex,aoindex));
      
      for(auto& gaussian:_gaussians){
        gaussian.contraction[contraction_index]/=norm;
      }
      aoindex+=numsubshell[i];
      contraction_index++;
    }
    
    return;
  }
  
  
  void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues, Eigen::Block< Eigen::MatrixX3d >& gradAOvalues, const tools::vec& grid_pos)const {

      // need position of shell
      const tools::vec center = grid_pos - _pos;
      const double center_x = center.getX();
      const double center_y = center.getY();
      const double center_z = center.getZ();
      const double distsq = center*center;



      // iterate over Gaussians in this shell
      for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

        const double alpha = itr->getDecay();
        const std::vector<double>& _contractions = itr->getContraction();

        const double _expofactor = itr->getPowfactor() * exp(-alpha * distsq);

        // split combined shells
        int _i_func = -1;
        int i_act;
        for (const char& single_shell:_type) {
          // single type shells
          if (single_shell == 'S') {
           
            i_act = _i_func + 1;
            AOvalues(i_act) += _contractions[0] * _expofactor; // s-function
            const double temp = _contractions[0] * -2.0 * alpha *_expofactor;
            gradAOvalues( i_act,0) += temp* center_x; // x gradient of s-function
            gradAOvalues( i_act,1) += temp* center_y; // y gradient of s-function
            gradAOvalues( i_act,2) += temp* center_z; // z gradient of s-function
            _i_func++;
          } else if (single_shell == 'P') {
            const double factor = 2. * sqrt(alpha) * _contractions[1];


            i_act = _i_func + 1; // Y 1,0
            double AOvalue = factor * center_z * _expofactor; // Y 1,0
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += -2. * center_x * alpha * AOvalue; // x gradient
            gradAOvalues( i_act,1) += -2. * center_y * alpha * AOvalue; // y gradient
            gradAOvalues( i_act,2) += factor * _expofactor - 2. * center_z * alpha * AOvalue; // z gradient

            i_act++; // Y 1,-1
            AOvalue = factor * center_y * _expofactor; // Y 1,-1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += -2. * center_x * alpha * AOvalue; // x gradient
            gradAOvalues( i_act,1) += factor * _expofactor - 2. * center_y * alpha * AOvalue; // y gradient
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue; // z gradient

            i_act++; // Y 1,1
            AOvalue = factor * center_x * _expofactor; // Y 1,1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor * _expofactor - 2. * center_x * alpha * AOvalue; // x gradient
            gradAOvalues( i_act,1) += -2. * center_y * alpha * AOvalue; // y gradient
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue; // z gradient

            _i_func += 3;
          } else if (single_shell == 'D') {
            const double factor = 2. * alpha * _contractions[2];
            const double factor_1 = factor / sqrt(3.);


            i_act = _i_func + 1; // Y 2,0
            double AOvalue = factor_1 * (3. * center_z * center_z - distsq) * _expofactor; // Y 2,0
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_1 * (-2. * center_x) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_1 * (-2. * center_y) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_1 * (4. * center_z) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 2,-1
            AOvalue = 2. * factor * (center_y * center_z) * _expofactor; // Y 2,-1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += -2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 2. * factor * center_z * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += 2. * factor * center_y * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 2,1
            AOvalue = 2. * factor * (center_x * center_z) * _expofactor; // Y 2,1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 2. * factor * center_z * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += -2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += 2. * factor * center_x * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 2,-2
            AOvalue = 2. * factor * (center_x * center_y) * _expofactor; // Y 2,-2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 2. * factor * center_y * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 2. * factor * center_x * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            i_act++; // Y 2,2
            AOvalue = factor * (center_x * center_x - center_y * center_y) * _expofactor; // Y 2,2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor * (2. * center_x) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor * (-2. * center_y) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            _i_func += 5;
          } else if (single_shell == 'F') {
            const double factor = 2. * pow(alpha, 1.5) * _contractions[3];
            const double factor_1 = factor * 2. / sqrt(15.);
            const double factor_2 = factor * sqrt(2.) / sqrt(5.);
            const double factor_3 = factor * sqrt(2.) / sqrt(3.);
            const double cx_cx = center_x*center_x;
            const double cx_cy = center_x*center_y;
            const double cx_cz = center_x*center_z;
            const double cy_cy = center_y*center_y;
            const double cy_cz = center_y*center_z;
            const double cz_cz = center_z*center_z;


            i_act = _i_func + 1; // Y 3,0
            double AOvalue = factor_1 * center_z * (5. * cz_cz - 3. * distsq) * _expofactor; // Y 3,0
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_1 * (-6. * cx_cz) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_1 * (-6. * cy_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            ///                      gradAOvalues( i_act,2) += factor_1 * (6. * cz_cz) * _expofactor - 2. * center_z * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_1 * 3. * (3. * cz_cz - distsq) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,-1
            AOvalue = factor_2 * center_y * (5. * cz_cz - distsq) * _expofactor; // Y 3,-1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_2 * (-2. * cx_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_2 * (4. * cz_cz - cx_cx - 3. * cy_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_2 * (8. * cy_cz) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,1
            AOvalue = factor_2 * center_x * (5. * cz_cz - distsq) * _expofactor; // Y 3,1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_2 * (4. * cz_cz - cy_cy - 3. * cx_cx) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_2 * (-2. * cx_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_2 * (8. * cx_cz) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,-2
            AOvalue = 4. * factor * center_x * center_y * center_z * _expofactor; // Y 3,-2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 4. * factor * cy_cz * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 4. * factor * cx_cz * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += 4. * factor * cx_cy * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,2
            AOvalue = 2. * factor * center_z * (cx_cx - cy_cy) * _expofactor; // Y 3,2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 2. * factor * (2. * cx_cz) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 2. * factor * (-2. * cy_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += 2. * factor * (cx_cx - cy_cy) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,-3
            AOvalue = factor_3 * center_y * (3. * cx_cx - cy_cy) * _expofactor; // Y 3,-3
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_3 * (6. * cx_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_3 * (3. * (cx_cx - cy_cy)) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            i_act++; // Y 3,3
            AOvalue = factor_3 * center_x * (cx_cx - 3. * cy_cy) * _expofactor; // Y 3,3
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_3 * (3. * (cx_cx - cy_cy)) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_3 * (-6. * cx_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            _i_func += 7;
          } else if (single_shell == 'G') {
            const double factor = 2. / sqrt(3.) * alpha * alpha * _contractions[4];
            const double factor_1 = factor / sqrt(35.);
            const double factor_2 = factor * 4. / sqrt(14.);
            const double factor_3 = factor * 2. / sqrt(7.);
            const double factor_4 = factor * 2. * sqrt(2.);
            const double cx_cx = center_x*center_x;
            const double cx_cy = center_x*center_y;
            const double cx_cz = center_x*center_z;
            const double cy_cy = center_y*center_y;
            const double cy_cz = center_y*center_z;
            const double cz_cz = center_z*center_z;


            i_act = _i_func + 1; // Y 4,0
            double AOvalue = factor_1 * (35. * cz_cz * cz_cz - 30. * cz_cz * distsq + 3. * distsq * distsq) * _expofactor; // Y 4,0
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_1 * 12. * center_x * (distsq - 5. * cz_cz) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_1 * 12. * center_y * (distsq - 5. * cz_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_1 * 16. * center_z * (5. * cz_cz - 3. * distsq) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,-1
            AOvalue = factor_2 * cy_cz * (7. * cz_cz - 3. * distsq) * _expofactor; // Y 4,-1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_2 * (-6. * center_x * cy_cz) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_2 * center_z * (4. * cz_cz - 3. * cx_cx - 9. * cy_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_2 * 3. * center_y * (5. * cz_cz - distsq) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,1
            AOvalue = factor_2 * cx_cz * (7. * cz_cz - 3. * distsq) * _expofactor; // Y 4,1
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_2 * center_z * (4. * cz_cz - 9. * cx_cx - 3. * cy_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_2 * (-6. * center_y * cx_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_2 * 3. * center_x * (5. * cz_cz - distsq) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,-2
            AOvalue = 2. * factor_3 * cx_cy * (7. * cz_cz - distsq) * _expofactor; // Y 4,-2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 2. * factor_3 * center_y * (6. * cz_cz - 3. * cx_cx - cy_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 2. * factor_3 * center_x * (6. * cz_cz - cx_cx - 3. * cy_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += 2. * factor_3 * 12. * center_z * cx_cy * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,2
            AOvalue = factor_3 * (cx_cx - cy_cy) * (7. * cz_cz - distsq) * _expofactor; // Y 4,2
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_3 * 4. * center_x * (3. * cz_cz - cx_cx) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_3 * 4. * center_y * (cy_cy - 3. * cz_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_3 * 12. * center_z * (cx_cx - cy_cy) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,-3
            AOvalue = factor_4 * cy_cz * (3. * cx_cx - cy_cy) * _expofactor; // Y 4,-3
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_4 * 6. * center_x * cy_cz * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_4 * 3. * center_z * (cx_cx - cy_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_4 * center_y * (3. * cx_cx - cy_cy) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,3
            AOvalue = factor_4 * cx_cz * (cx_cx - 3. * cy_cy) * _expofactor; // Y 4,3
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor_4 * 3. * center_z * (cx_cx - cy_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor_4 * (-6. * center_y * cx_cz) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += factor_4 * center_x * (cx_cx - 3. * cy_cy) * _expofactor - 2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,-4
            AOvalue = 4. * factor * cx_cy * (cx_cx - cy_cy) * _expofactor; // Y 4,-4
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += 4. * factor * center_y * (3. * cx_cx - cy_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += 4. * factor * center_x * (cx_cx - 3. * cy_cy) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            i_act++; // Y 4,4
            AOvalue = factor * (cx_cx * cx_cx - 6. * cx_cx * cy_cy + cy_cy * cy_cy) * _expofactor; // Y 4,4
            AOvalues(i_act) += AOvalue;
            gradAOvalues( i_act,0) += factor * 4. * center_x * (cx_cx - 3. * cy_cy) * _expofactor - 2. * center_x * alpha * AOvalue;
            gradAOvalues( i_act,1) += factor * 4. * center_y * (cy_cy - 3. * cx_cx) * _expofactor - 2. * center_y * alpha * AOvalue;
            gradAOvalues( i_act,2) += -2. * center_z * alpha * AOvalue;

            _i_func += 9;
            ;
          } else if (single_shell == 'H') {
            std::cerr << "H functions not implemented in AOeval at the moment!" << std::endl;
            exit(1);
          } else {
            std::cerr << "Single shell type" << single_shell << " not known " << std::endl;
            exit(1);
          }
        }
      } // contractions
      return;
    }

           
     
void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues, const tools::vec& grid_pos )const{

           
            // need position of shell
             const tools::vec center=grid_pos-this->_pos;
             const double center_x = center.getX();
             const double center_y = center.getY();
             const double center_z = center.getZ();
             const double distsq =  center*center;

            
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double alpha = itr->getDecay();
                const std::vector<double>& _contractions = itr->getContraction();

                double _expofactor =itr->getPowfactor() * exp(-alpha * distsq);

                // split combined shells
                int _i_func = -1;

                for (const char& single_shell :_type) {
                    // single type shells
                    if (single_shell == 'S') {
                        AOvalues( _i_func + 1) += _contractions[0] * _expofactor; // s-function        
                        _i_func++;
                    }
                    else if (single_shell == 'P') {
                        double factor = 2.*sqrt(alpha)*_contractions[1];
                      AOvalues( _i_func + 1) += factor * center_z* _expofactor; // Y 1,0
                      AOvalues( _i_func + 2) += factor * center_y* _expofactor; // Y 1,-1
                      AOvalues(_i_func + 3) += factor * center_x* _expofactor; // Y 1,1
                      _i_func += 3;
                    }
                    else if (single_shell == 'D') {
                       double factor = 2.*alpha*_contractions[2];
                      double factor_1 =  factor/sqrt(3.);
                      AOvalues( _i_func + 1) += factor_1 * (3.*center_z*center_z - distsq) * _expofactor; // Y 2,0
                      AOvalues( _i_func + 2) += 2.*factor * (center_y*center_z) * _expofactor; // Y 2,-1
                      AOvalues( _i_func + 3) += 2.*factor * (center_x*center_z) * _expofactor; // Y 2,1
                      AOvalues( _i_func + 4) += 2.*factor * (center_x*center_y) * _expofactor; // Y 2,-2
                      AOvalues(_i_func + 5) += factor * (center_x*center_x - center_y*center_y) * _expofactor; // Y 2,2
                      _i_func += 5;
                    }
                    else if (single_shell == 'F') {
                      double factor = 2.*pow(alpha,1.5)*_contractions[3];
                      double factor_1 = factor*2./sqrt(15.);
                      double factor_2 = factor*sqrt(2.)/sqrt(5.);
                      double factor_3 = factor*sqrt(2.)/sqrt(3.);
                      double cx_cx = center_x*center_x;
                      double cy_cy = center_y*center_y;
                      double cz_cz = center_z*center_z;

                      AOvalues(_i_func + 1) += factor_1 * center_z * (5.*cz_cz - 3.*distsq) * _expofactor; // Y 3,0
                      AOvalues(_i_func + 2) += factor_2 * center_y * (5.*cz_cz - distsq) * _expofactor; // Y 3,-1
                      AOvalues( _i_func + 3) += factor_2 * center_x * (5.*cz_cz - distsq) * _expofactor; // Y 3,1
                      AOvalues( _i_func + 4) += 4.*factor * center_x * center_y * center_z * _expofactor; // Y 3,-2
                      AOvalues( _i_func + 5) += 2.*factor * center_z * (cx_cx - cy_cy) * _expofactor; // Y 3,2
                      AOvalues( _i_func + 6) += factor_3 * center_y * (3.*cx_cx - cy_cy) * _expofactor; // Y 3,-3
                      AOvalues( _i_func + 7) += factor_3 * center_x * (cx_cx - 3.*cy_cy) * _expofactor; // Y 3,3

                      _i_func += 7;
                    }
                    else if (single_shell == 'G') {
                      double factor = 2./sqrt(3.)*alpha*alpha*_contractions[4];
                      double factor_1 = factor/sqrt(35.);
                      double factor_2 = factor*4./sqrt(14.);
                      double factor_3 = factor*2./sqrt(7.);
                      double factor_4 = factor*2.*sqrt(2.);
                      double cx_cx = center_x*center_x;
                      double cx_cy = center_x*center_y;
                      double cx_cz = center_x*center_z;
                      double cy_cy = center_y*center_y;
                      double cy_cz = center_y*center_z;
                      double cz_cz = center_z*center_z;

                      AOvalues( _i_func + 1) += factor_1 * (35.*cz_cz*cz_cz - 30.*cz_cz*distsq + 3.*distsq*distsq) * _expofactor; // Y 4,0
                      AOvalues(_i_func + 2) += factor_2 * cy_cz * (7.*cz_cz - 3.*distsq) * _expofactor; // Y 4,-1
                      AOvalues( _i_func + 3) += factor_2 * cx_cz * (7.*cz_cz - 3.*distsq) * _expofactor; // Y 4,1
                      AOvalues( _i_func + 4) += 2.*factor_3 * cx_cy * (7.*cz_cz - distsq) * _expofactor; // Y 4,-2
                      AOvalues(_i_func + 5) += factor_3 * (cx_cx - cy_cy) * (7.*cz_cz - distsq) * _expofactor; // Y 4,2
                      AOvalues( _i_func + 6) += factor_4 * cy_cz * (3.*cx_cx - cy_cy) * _expofactor; // Y 4,-3
                      AOvalues( _i_func + 7) += factor_4 * cx_cz * (cx_cx - 3.*cy_cy) * _expofactor; // Y 4,3
                      AOvalues( _i_func + 8) += 4.*factor * cx_cy * (cx_cx - cy_cy) * _expofactor; // Y 4,-4
                      AOvalues(_i_func + 9) += factor * (cx_cx*cx_cx - 6.*cx_cx*cy_cy + cy_cy*cy_cy) * _expofactor; // Y 4,4

                      _i_func += 9;
                    }
                    else{
                        std::cerr << "Single shell type"<<single_shell<<" not known " << std::endl;
                        exit(1);
                    }
                }
            } // contractions
            return;
        }
        



}}
