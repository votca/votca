/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


namespace votca { namespace xtp {
    
 

int AOShell::detlmax( string shell_type ) {
    int _lmax;
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _lmax = 0;}
       if ( shell_type == "P" ){ _lmax = 1;}
       if ( shell_type == "D" ){ _lmax = 2;}
       if ( shell_type == "F" ){ _lmax = 3;}
       if ( shell_type == "G" ){ _lmax = 4;}
    } else {
        // for combined shells check all contributions
        _lmax = 0;
        for( unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell =    string( shell_type, i, 1 );
            int _test = this->detlmax( local_shell  );
            if ( _test > _lmax ) { _lmax = _test;} 
        }
    }

    return _lmax;
}


    
void AOShell::EvalAOGradspace(ub::matrix_range<ub::matrix<double> >& gradAOvalues, double x, double y, double z, string type ){

     // need type of shell
     string shell_type = this->_type;
     // need position of shell
     double center_x = x - this->_pos.getX();
     double center_y = y - this->_pos.getY();
     double center_z = z - this->_pos.getZ();
     double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
     const double pi = boost::math::constants::pi<double>();
            
            
  typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
           // std::vector< std::vector<double> > _AOGradevaluated(3);
         
            // split combined shells
            int _i_func = -1;
            int i_act;
            for (unsigned i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[0] * -2.0*alpha*center_x*_expofactor; // x gradient of s-function
                    gradAOvalues(i_act,1) += _contractions[0] * -2.0*alpha*center_y*_expofactor; // y gradient of s-function
                    gradAOvalues(i_act,2) += _contractions[0] * -2.0*alpha*center_z*_expofactor; // z gradient of s-function
                    _i_func = i_act;
                }
                if ( single_shell == "P") {
                    
                    // px -functions
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_x*center_x)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // z gradient 
                    
                    // py -functions
                    i_act = _i_func +2;
                    gradAOvalues(i_act,0) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_y*center_y) *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z*_expofactor; // z gradient 
                    
                     // pz -functions
                    i_act = _i_func +3;
                    gradAOvalues(i_act,0) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_z*center_z)*_expofactor; // z gradient                    
                    _i_func = i_act;

                }
                if ( single_shell == "D") {
                             
                    // dxz function
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * (center_z - 2.0*alpha*center_x*center_x*center_z)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * -8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha * (center_x - 2.0*alpha*center_x*center_z*center_z)*_expofactor; // z gradient 

                    // dyz function
                    i_act = _i_func+2;                    
                    gradAOvalues(i_act,0) += _contractions[2] * -8.0*alpha*alpha * center_x*center_y*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * 4.0*alpha * (center_z - 2.0*alpha*center_y*center_y*center_z)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha * (center_y - 2.0*alpha*center_y*center_z*center_z)*_expofactor; // z gradient 

                    // dxy function
                    i_act = _i_func+3;                    
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * (center_y - 2.0*alpha*center_x*center_x*center_y)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * 4.0*alpha * (center_x - 2.0*alpha*center_x*center_y*center_y)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) +=_contractions[2] * -8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // z gradient 

                    // d3z2-r2-function
                    i_act = _i_func+4;                    
                    gradAOvalues(i_act,0) +=_contractions[2] * -4.0*alpha/sqrt(3.0) * center_x * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,1) +=_contractions[2] * -4.0*alpha/sqrt(3.0) * center_y * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha/sqrt(3.0) * center_z * (2.0 - alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    
                    // dx2-y2-function
                    i_act = _i_func+5;                    
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * center_x * (1.0 - alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,1) +=_contractions[2] * -4.0*alpha * center_y * (1.0 + alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,2) += _contractions[2] * -4.0*alpha*alpha * center_z * (center_x*center_x - center_y*center_y) * _expofactor;


                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
            }
            }// contractions


        }
    /*
       inline void AOShell::EvalAOGradspace(ub::matrix<double>& gradAOvalues, double x, double y, double z, string type ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            double alpha = this->_gaussians[0]->decay; // only uncontracted for testing
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
            std::vector< std::vector<double> > _AOGradevaluated(3);
         
            // split combined shells
            int _i_func = -1;
            int i_act;
            for (int i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) =-2.0*alpha*center_x*_expofactor; // x gradient of s-function
                    gradAOvalues(i_act,1) =-2.0*alpha*center_y*_expofactor; // y gradient of s-function
                    gradAOvalues(i_act,2) =-2.0*alpha*center_z*_expofactor; // z gradient of s-function
                    _i_func = i_act;
                }
                if ( single_shell == "P") {
                    
                    // px -functions
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_x*center_x)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // z gradient 
                    
                    // py -functions
                    i_act = _i_func +2;
                    gradAOvalues(i_act,0) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_y*center_y) *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z*_expofactor; // z gradient 
                    
                     // pz -functions
                    i_act = _i_func +3;
                    gradAOvalues(i_act,0) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_z*center_z)*_expofactor; // z gradient                    
                    _i_func = i_act;

                }
                if ( single_shell == "D") {
                             // dxz, dyz, dxy, d3z2-r2, dx2-y2
                    // dxz function
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) = 4.0*alpha * (center_z - 2.0*alpha*center_x*center_x*center_z)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) = 4.0*alpha * (center_x - 2.0*alpha*center_x*center_z*center_z)*_expofactor; // z gradient 

                    // dyz function
                    i_act = _i_func+2;                    
                    gradAOvalues(i_act,0) =-8.0*alpha*alpha * center_x*center_y*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) = 4.0*alpha * (center_z - 2.0*alpha*center_y*center_y*center_z)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) = 4.0*alpha * (center_y - 2.0*alpha*center_y*center_z*center_z)*_expofactor; // z gradient 

                    // dxy function
                    i_act = _i_func+3;                    
                    gradAOvalues(i_act,0) = 4.0*alpha * (center_y - 2.0*alpha*center_x*center_x*center_y)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) = 4.0*alpha * (center_x - 2.0*alpha*center_x*center_y*center_y)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // z gradient 

                    // d3z2-r2-function
                    i_act = _i_func+4;                    
                    gradAOvalues(i_act,0) =-4.0*alpha/sqrt(3.0) * center_x * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,1) =-4.0*alpha/sqrt(3.0) * center_y * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,2) = 4.0*alpha/sqrt(3.0) * center_z * (2.0 - alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    
                    // dx2-y2-function
                    i_act = _i_func+5;                    
                    gradAOvalues(i_act,0) = 4.0*alpha * center_x * (1.0 - alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,1) =-4.0*alpha * center_y * (1.0 + alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,2) =-4.0*alpha*alpha * center_z * (center_x*center_x - center_y*center_y) * _expofactor;


                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
            }



        } */
    
        
       
void AOShell::EvalAOIntegral(ub::matrix_range<ub::matrix<double> >& AOvalues){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
      //      double center_x = x - this->_pos.getX();
    //        double center_y = y - this->_pos.getY();
  //          double center_z = z - this->_pos.getZ();
            // need decay constant
            
//            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            
            
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _factor = pow(2.0 * pi / alpha, 0.75) ;

                // split combined shells
                int _i_func = -1;
                //int i_act;
                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _factor; // s-function
                        _i_func++;
                    }
                }
            } // contractions

        }
       
       
       
void AOShell::EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, ub::matrix_range<ub::matrix<double> >& gradAOvalues, double x, double y, double z ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            
            
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _expofactor = pow(2.0 * alpha / pi, 0.75) * exp(-alpha * distsq);

                // split combined shells
                int _i_func = -1;
                int i_act;
                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _expofactor; // s-function

                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[0] * -2.0 * alpha * center_x*_expofactor; // x gradient of s-function
                            gradAOvalues(1, i_act) += _contractions[0] * -2.0 * alpha * center_y*_expofactor; // y gradient of s-function
                            gradAOvalues(2, i_act) += _contractions[0] * -2.0 * alpha * center_z*_expofactor; // z gradient of s-function
                        
                        _i_func++;
                    }
                    if (single_shell == "P") {
                        AOvalues(0, _i_func + 1) += _contractions[1] * 2.0 * sqrt(alpha) * center_x*_expofactor; // px-function
                        AOvalues(0, _i_func + 2) += _contractions[1] * 2.0 * sqrt(alpha) * center_y*_expofactor; // py-function
                        AOvalues(0, _i_func + 3) += _contractions[1] * 2.0 * sqrt(alpha) * center_z*_expofactor; // pz-function
                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_x * center_x) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_y*_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_z*_expofactor; // z gradient 

                            // py -functions
                            i_act = _i_func + 2;
                            gradAOvalues(0, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_y*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_y * center_y) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_y * center_z*_expofactor; // z gradient 

                            // pz -functions
                            i_act = _i_func + 3;
                            gradAOvalues(0, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_z*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_y * center_z *_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_z * center_z) * _expofactor; // z gradient      



                        _i_func += 3;
                    }
                    if (single_shell == "D") {
                        // dxz, dyz, dxy, d3z2-r2, dx2-y2
                        AOvalues(0, _i_func + 1) += _contractions[2] * 4.0 * alpha * center_x * center_z*_expofactor; // dxz-function
                        AOvalues(0, _i_func + 2) += _contractions[2] * 4.0 * alpha * center_y * center_z*_expofactor; // dyz-function
                        AOvalues(0, _i_func + 3) += _contractions[2] * 4.0 * alpha * center_x * center_y*_expofactor; // dxy-function
                        AOvalues(0, _i_func + 4) += _contractions[2] * 2.0 * alpha / sqrt(3.0)*(3.0 * center_z * center_z - distsq) * _expofactor; // d3z2-r2-function
                        AOvalues(0, _i_func + 5) += _contractions[2] * 2.0 * alpha * (center_x * center_x - center_y * center_y) * _expofactor; // dx2-y2-function


                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * (center_z - 2.0 * alpha * center_x * center_x * center_z) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha * (center_x - 2.0 * alpha * center_x * center_z * center_z) * _expofactor; // z gradient 

                            // dyz function
                            i_act = _i_func + 2;
                            gradAOvalues(0, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * 4.0 * alpha * (center_z - 2.0 * alpha * center_y * center_y * center_z) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha * (center_y - 2.0 * alpha * center_y * center_z * center_z) * _expofactor; // z gradient 

                            // dxy function
                            i_act = _i_func + 3;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * (center_y - 2.0 * alpha * center_x * center_x * center_y) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * 4.0 * alpha * (center_x - 2.0 * alpha * center_x * center_y * center_y) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // z gradient 

                            // d3z2-r2-function
                            i_act = _i_func + 4;
                            gradAOvalues(0, i_act) += _contractions[2] * -4.0 * alpha / sqrt(3.0) * center_x * (1.0 + alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;
                            gradAOvalues(1, i_act) += _contractions[2] * -4.0 * alpha / sqrt(3.0) * center_y * (1.0 + alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha / sqrt(3.0) * center_z * (2.0 - alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;

                            // dx2-y2-function
                            i_act = _i_func + 5;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * center_x * (1.0 - alpha * (center_x * center_x - center_y * center_y)) * _expofactor;
                            gradAOvalues(1, i_act) += _contractions[2] * -4.0 * alpha * center_y * (1.0 + alpha * (center_x * center_x - center_y * center_y)) * _expofactor;
                            gradAOvalues(2, i_act) += _contractions[2] * -4.0 * alpha * alpha * center_z * (center_x * center_x - center_y * center_y) * _expofactor;



                        
                        
                        
                        _i_func += 5;
                    }
                    if (single_shell == "F") {
                        cerr << " F functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                    if (single_shell == "G") {
                        cerr << " G functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                }
            } // contractions

        }
           
           
           
void AOShell::EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, double x, double y, double z ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();

            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _expofactor = pow(2.0 * alpha / pi, 0.75) * exp(-alpha * distsq);

                // split combined shells
                int _i_func = -1;

                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _expofactor; // s-function        
                        _i_func++;
                    }
                    if (single_shell == "P") {
                        AOvalues(0, _i_func + 1) += _contractions[1] * 2.0 * sqrt(alpha) * center_x*_expofactor; // px-function
                        AOvalues(0, _i_func + 2) += _contractions[1] * 2.0 * sqrt(alpha) * center_y*_expofactor; // py-function
                        AOvalues(0, _i_func + 3) += _contractions[1] * 2.0 * sqrt(alpha) * center_z*_expofactor; // pz-function

                        _i_func += 3;
                    }
                    if (single_shell == "D") {
                        // dxz, dyz, dxy, d3z2-r2, dx2-y2
                        AOvalues(0, _i_func + 1) += _contractions[2] * 4.0 * alpha * center_x * center_z*_expofactor; // dxz-function
                        AOvalues(0, _i_func + 2) += _contractions[2] * 4.0 * alpha * center_y * center_z*_expofactor; // dyz-function
                        AOvalues(0, _i_func + 3) += _contractions[2] * 4.0 * alpha * center_x * center_y*_expofactor; // dxy-function
                        AOvalues(0, _i_func + 4) += _contractions[2] * 2.0 * alpha / sqrt(3.0)*(3.0 * center_z * center_z - distsq) * _expofactor; // d3z2-r2-function
                        AOvalues(0, _i_func + 5) += _contractions[2] * 2.0 * alpha * (center_x * center_x - center_y * center_y) * _expofactor; // dx2-y2-function                             
                        _i_func += 5;
                    }
                    if (single_shell == "F") {
                        cerr << " F functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                    if (single_shell == "G") {
                        cerr << " G functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                }
            } // contractions

        }
        
           
           
           
        
           
        
               
  /*         inline void AOShell::EvalAOspace(ub::matrix<double>& AOvalues, double x, double y, double z, string type ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            double alpha = this->_gaussians[0]->decay; // only uncontracted for testing
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
            // split combined shells
            int _i_func = -1;
            for (int i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    AOvalues(_i_func +1 ,0) = _expofactor; // s-function
                    _i_func++;
                }
                if ( single_shell == "P") {
                    AOvalues(_i_func +1,0) = 2.0*sqrt(alpha)*center_x*_expofactor; // px-function
                    AOvalues(_i_func +2,0) = 2.0*sqrt(alpha)*center_y*_expofactor; // py-function
                    AOvalues(_i_func +3,0) = 2.0*sqrt(alpha)*center_z*_expofactor; // pz-function
                    _i_func += 3;
                }
                if ( single_shell == "D") {
                    // dxz, dyz, dxy, d3z2-r2, dx2-y2
                    AOvalues(_i_func +1,0) = 4.0*alpha*center_x*center_z*_expofactor; // dxz-function
                    AOvalues(_i_func +2,0) = 4.0*alpha*center_y*center_z*_expofactor; // dyz-function
                    AOvalues(_i_func +3,0) = 4.0*alpha*center_x*center_y*_expofactor; // dxy-function
                    AOvalues(_i_func +4,0) = 2.0*alpha/sqrt(3.0)*(3.0*center_z*center_z - distsq)*_expofactor; // d3z2-r2-function
                    AOvalues(_i_func +5,0) = 2.0*alpha*(center_x*center_x - center_y*center_y)*_expofactor; // dx2-y2-function
                    _i_func += 5;
                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOeval at the moment!" << endl;
                    exit(1);
                }
            }   


        } */
        
        


}}
