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

#ifndef __XTP_EXCHANGE_CORRELATION__H
#define	__XTP_EXCHANGE_CORRELATION__H


#include <boost/math/constants/constants.hpp>
using namespace std;


namespace votca {
    namespace xtp {

        class ExchangeCorrelation {
        public:

            ExchangeCorrelation() {
            };

            void getXC(int type, const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& f, double& df_drho, double& df_dsigma);



        private:

            void set_rs(const double& rho, double& rs, double& drs_dr);
            void set_fzeta(const double& zeta, double& fzeta, double& dfzeta_dzeta);
            void set_s(const double& density, double& gradNorm, double& s, double& ds_drho, double& ds_dnrho);
            void setNormGradient(const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& norm, double& drho_dX_normal, double& drho_dY_normal, double& drho_dZ_normal);
            void set_gphi(const double& zeta, double& gphi, double& dgphi_dzeta);
            void set_t(const double& rho, double& gphi, double& nrho, double& t, double& dt_drho, double& dt_dgphi, double& dt_dnrho  );

            void evalLDA(const double& rho, double& f, double& df_drho);
            void evalLDA_exchange(const double& rho, double& fx, double& dfx_drho, const double& ex_par);
            void evalLDA_correlation(const double& rho, double& fc, double& dfc_drho);

            void evalPW91_exchange(const double& rho, double& fx, double& dfx_drho);
            void evalPW91_correlation_G1( const double& rs, double& G, double& dG_drs, const double& p);
            void evalPW91_correlation_G2( const double& rs, double& G, double& dG_drs, const double& p);
            void evalPW91_correlation_G3( const double& rs, double& G, double& dG_drs, const double& p);
            void evalPW91_correlation_G( const double& rs, double& G, double& dG_drs, const double& A, const double& alpha1, const double& beta1, const double& beta2, const double& beta3, const double& beta4, const double& p);
            void evalPW91_correlation_int( const double& rho, const double& zeta, double& ec, double& dec_drho, double& dec_dzeta );

            //void evalPBE(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& f, double& df_drho, double& df_drho_dX, double& df_drho_dY, double& df_drho_dZ);
            void evalPBE(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& f, double& df_drho, double& df_dsigma);
            void evalPBE_exchange(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& fx, double& dfx_drho, double& dfx_drho_dX, double& dfx_drho_dY, double& dfx_drho_dZ);
            void evalPBE_correlation(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& fc, double& dfc_drho, double& dfc_drho_dX,double& dfc_drho_dY,double& dfc_drho_dZ);
            void evalPBE_correlation_H(const double& rho,  double& zeta,  double& nrho,  double& eps,  double& deps_drho,  double& deps_dzeta, double& H, double& dH_drho, double& dH_dzeta, double& dH_dnrho);

        };


        // Main interface function to calculate XC potential and it's gradient, depending on local rho and gradient of rho        

        inline void ExchangeCorrelation::getXC(int type, const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& f, double& df_drho, double& df_dsigma) {

            if (type == -1) evalLDA(rho, f, df_drho);
            else if (type == -2) evalPBE(rho, drho_dX, drho_dY, drho_dZ, f, df_drho, df_dsigma );


        }

        /* ======= GGA PBE ======== */

        inline void ExchangeCorrelation::evalPBE(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& f, double& df_drho, double& df_dsigma) {

            double fx = 0.0;
            double dfx_drho = 0.0;
            double dfx_drho_dX = 0.0;
            double dfx_drho_dY = 0.0;
            double dfx_drho_dZ = 0.0;
            evalPBE_exchange(rho, drho_dX, drho_dY, drho_dZ, fx, dfx_drho, dfx_drho_dX, dfx_drho_dY, dfx_drho_dZ);

	    // cout  <<  " VOTCA exchange vrho " << dfx_drho << " v_grad_x " << dfx_drho_dX << " v_grad_y " << dfx_drho_dY  << " v_grad_z " << dfx_drho_dZ << endl;
            // cout  <<  " VOTCA exchange vrho " << dfx_drho << " v_sigma " << dfx_drho_dX/drho_dX/2 << endl;
            
            double fc = 0.0;
            double dfc_drho = 0.0;
            double dfc_drho_dX = 0.0;
            double dfc_drho_dY = 0.0;
            double dfc_drho_dZ = 0.0;
            evalPBE_correlation(rho, drho_dX, drho_dY, drho_dZ, fc, dfc_drho, dfc_drho_dX, dfc_drho_dY, dfc_drho_dZ);

	    //cout  <<  " VOTCA correlation vrho " << dfc_drho << " v_grad_x " << dfc_drho_dX << " v_grad_y " << dfc_drho_dY  << " v_grad_z " << dfc_drho_dZ << endl ;
            // cout  <<  " VOTCA correlation vrho " << dfc_drho << " v_sigma_x " << dfc_drho_dX/drho_dX/2 << endl;
            
            f          = fx          + fc;
            df_drho    = dfx_drho    + dfc_drho;
            df_dsigma  = 0.5/drho_dX * ( dfx_drho_dX + dfc_drho_dX ) ;
            //df_drho_dX = dfx_drho_dX + dfc_drho_dX;
            //df_drho_dY = dfx_drho_dY + dfc_drho_dY;
            //df_drho_dZ = dfx_drho_dZ + dfc_drho_dZ;
            
            
        }




        // GGA PBE EXCHANGE

        inline void ExchangeCorrelation::evalPBE_exchange(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& fx, double& dfx_drho, double& dfx_drho_dX, double& dfx_drho_dY, double& dfx_drho_dZ) {

            double norm = 0.0;
            double drho_dX_normal = 0.0;
            double drho_dY_normal = 0.0;
            double drho_dZ_normal = 0.0;


            if (rho <= 0.0) {
                return ;
            } else {
                // exchange part
                // constants
                const double mu = 0.21951;
                const double kappa = 0.804;

                // get PW91 exchange
                double fx_LDA = 0.0;
                double dfx_drho_LDA = 0.0;
                evalPW91_exchange(rho, fx_LDA, dfx_drho_LDA);

                double s;
                double ds_dr;
                double ds_dnrho;
                // prep norm of gradient and gradient normal vector
                setNormGradient(drho_dX, drho_dY, drho_dZ, norm, drho_dX_normal, drho_dY_normal, drho_dZ_normal);
                set_s(rho, norm, s, ds_dr, ds_dnrho);
                
                
                if (std::abs(s) < 1.e-10) {
                    
                    fx = fx_LDA;
                    dfx_drho = dfx_drho_LDA;

                } else {

                    double F_denom = 1.0 + mu * s * s / kappa;
                    double F = 1.0 + kappa - kappa / F_denom;
                    double dF_ds = 2.0 * mu * s / (F_denom * F_denom);
                    
                   // cout << "F_denom " << F_denom << " F " << F << " dF_ds " << dF_ds << endl;
                    
                    fx          = F * fx_LDA;
                    dfx_drho    = dfx_drho_LDA * F + fx_LDA * dF_ds * ds_dr;
                    dfx_drho_dX = fx_LDA * dF_ds * ds_dnrho * drho_dX_normal;
                    dfx_drho_dY = fx_LDA * dF_ds * ds_dnrho * drho_dY_normal;
                    dfx_drho_dZ = fx_LDA * dF_ds * ds_dnrho * drho_dZ_normal;
      
                    
                }

            } // exchange part

        } // evalPBE_exchange 


        // GGA PBE CORRELATION
        inline void ExchangeCorrelation::evalPBE_correlation(const double& rho, const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& fc, double& dfc_drho, double& dfc_drho_dX, double& dfc_drho_dY, double& dfc_drho_dZ) {

            double norm = 0.0;
            double drho_dX_normal = 0.0;
            double drho_dY_normal = 0.0;
            double drho_dZ_normal = 0.0;

            // prep norm of gradient and gradient normal vector
            setNormGradient(drho_dX, drho_dY, drho_dZ, norm, drho_dX_normal, drho_dY_normal, drho_dZ_normal);
            
            // again get something from PW91
            double zeta = 0.0;
            double eps = 0.0;
            double deps_drho = 0.0;
            double deps_dzeta = 0.0;
            evalPW91_correlation_int(rho,zeta,eps,deps_drho,deps_dzeta);
            
            // and H
            double H        = 0.0;
            double dH_drho  = 0.0;
            double dH_dzeta = 0.0;
            double dH_dnrho = 0.0;
            evalPBE_correlation_H(rho,zeta,  norm,  eps,  deps_drho,  deps_dzeta,H,dH_drho,dH_dzeta,dH_dnrho);

            //cout << "H " << H << " dh_drho " << dH_drho << " dH_dzeta  " << dH_dzeta << " dH_dnrho " << dH_dnrho << endl;
            
            
            // finallY
            fc       = rho * (eps + H);
            dfc_drho = (eps + H) + rho * (deps_drho + dH_drho);
            
            //cout << "rho " << rho << " eps " << eps << " deps_drho  " << deps_drho << " dH_drho " << dH_drho << endl;
            
            dfc_drho_dX = rho * dH_dnrho * drho_dX_normal;
            dfc_drho_dY = rho * dH_dnrho * drho_dY_normal;
            dfc_drho_dZ = rho * dH_dnrho * drho_dZ_normal;


            
        } // evalPBE_correlation
        
        
        
        inline void ExchangeCorrelation::evalPBE_correlation_H(const double& rho,  double& zeta,  double& nrho,  double& eps, double& deps_drho, double& deps_dzeta, double& H, double& dH_drho, double& dH_dzeta, double& dH_dnrho){
            
        
            const double beta = 0.066725;
            const double gamma = 0.031091;
            
            double gphi;
            double dgphi_dzeta;
            set_gphi(zeta,gphi,dgphi_dzeta);
          
            
            double t;
            double dt_drho;
            double dt_dgphi;
            double dt_dnrho;
            set_t(rho, gphi, nrho,  t,  dt_drho,  dt_dgphi, dt_dnrho  );


            if (std::abs(t) < 1e-10) {
                H = 0.0;
                dH_drho = 0.0;
                dH_dzeta = 0.0;
                dH_dnrho = 0.0;
                return;
            }
            
            // A part
            
            
            double Aea = - eps / (gamma * pow(gphi,3) ); // there was some Hartree in here, but we didn't convert eps to Hartree before!
            double Aexp = exp(Aea);
            double Adenom;
            if ( std::abs(Aea) < sqrt( std::numeric_limits<double>::epsilon() ) ){
                Adenom = Aea;
            } else {
                Adenom = Aexp - 1.0;
            }

            double A = beta / gamma / Adenom;
            double dAea_deps = Aea / eps;
            double dAea_dgphi = -3.0 * Aea / gphi;
            double dA_dAea = -beta / gamma / pow(Adenom, 2) * Aexp;
            double dA_deps = dA_dAea * dAea_deps;
            double dA_dgphi = dA_dAea * dAea_dgphi;

            // --- Hb

            double Hb = A * pow(t, 2);
            double dHb_dA = pow(t, 2);
            double dHb_dt = 2.0 * A*t;

            // --- Ha

            double Ha_denom = 1.0 + Hb + Hb*Hb;
            double Ha = beta / gamma * pow(t, 2) * (1.0 + Hb) / Ha_denom;

            double dHa_dHb = -beta / gamma * pow(t, 2) * Hb * (2.0 + Hb) / (Ha_denom * Ha_denom);
            double dHa_dt_noHb = 2.0 * Ha / t;

            double dHa_dt = dHa_dt_noHb + dHa_dHb * dHb_dt;
            double dHa_dA = dHa_dHb * dHb_dA;

            // --- H

            H = gamma * pow(gphi, 3) * log(1.0 + Ha); // there was a Hartree here

            double dH_dHa = gamma * pow(gphi, 3) / (1.0 + Ha); // there was a Hartree here

            double dH_dA = dH_dHa * dHa_dA;
            double dH_dt = dH_dHa * dHa_dt;

            double dH_dgphi = 3.0 * H / gphi + dH_dA * dA_dgphi + dH_dt*dt_dgphi;
            double dH_deps = dH_dA * dA_deps;

            dH_dzeta = dH_dgphi * dgphi_dzeta + dH_deps * deps_dzeta;
            dH_drho = dH_deps * deps_drho + dH_dt * dt_drho;
            
            // cout << " rho " << rho << " dH_deps " << dH_deps << " deps_drho " << deps_drho << " dH_dt " << dH_dt << " dt_drho " << dt_drho << endl; 
            
            dH_dnrho = dH_dt * dt_dnrho;
            
            
            
            
        }
        
        
        
       inline void ExchangeCorrelation::evalPW91_correlation_int(const double& rho, const double& zeta, double& ec, double& dec_drho, double& dec_dzeta){
            
                 const double f0_second = 1.709920934161365447;

                 // get zeta derived values ( all 0 for PBE )
                 double fzeta = 0.0 ;
                 double dfzeta_dzeta = 0.0 ;
                 if ( zeta != 0.0 ) set_fzeta(zeta,fzeta,dfzeta_dzeta);
                 
                 // get rs
                 double rs;
                 double drs_drho;
                 set_rs(rho, rs, drs_drho);
                 
                 // first evaluation of PW91 G function 
                 double ec0;
                 double dec0_drs;
                 double p = 1.0;
                 evalPW91_correlation_G1(rs,ec0,dec0_drs,p);
                 double dec0_drho = dec0_drs * drs_drho;

                 if ( zeta != 0.0 ){
                 
                 // second evaluation of PW91 G function
                 double ec1;
                 double dec1_drs;
                 evalPW91_correlation_G2(rs,ec1,dec1_drs,p);
                 double dec1_drho = dec1_drs * drs_drho;
                 
                 // third evaluation of PW91 G function
                 double eca;
                 double deca_drs;
                 evalPW91_correlation_G3(rs,eca,deca_drs,p);
                 eca = - eca; 
                 deca_drs = - deca_drs;
                 double deca_drho = deca_drs * drs_drho;
                 
                 // now sum all contributions (weighted)
                 double afac = fzeta / f0_second * (1.0 - pow(zeta,4));
                 double dafac_dzeta = dfzeta_dzeta / f0_second * (1.0 - pow(zeta,4) ) - fzeta / f0_second * 4.0*pow(zeta,3);
                 double pfac = fzeta * pow(zeta,4);
                 double dpfac_dzeta = dfzeta_dzeta * pow(zeta,4) + fzeta * 4.0*pow(zeta,3);

                 ec = ec0 + eca*afac + (ec1 - ec0)*pfac;
                 dec_drho = dec0_drho + deca_drho*afac + (dec1_drho - dec0_drho)*pfac;
                 dec_dzeta = eca*dafac_dzeta + (ec1 - ec0)*dpfac_dzeta;

                 } else {
                     
                     ec = ec0;
                     dec_drho = dec0_drho;
                     dec_dzeta = 0.0;
                     
                     
                 }
                 

        } 



       inline void ExchangeCorrelation::evalPW91_correlation_G(const double& rs, double& G, double& dG_drs, const double& A, const double& alpha1, const double& beta1, const double& beta2, const double& beta3, const double& beta4, const double& p){


             double rs_1 = sqrt(rs);
             double rs_3 = rs * rs_1;
             double rs_p = pow(rs, p);
             double rs_p1 = rs * rs_p;

             double ln_denom = 2.0 * A * (beta1 * rs_1 + beta2 * rs + beta3 * rs_3 + beta4 * rs_p1);
             double ln_arg = 1.0 + 1.0 / ln_denom;
             double ln_res = log(ln_arg);
             double pre_fac = -2.0 * A * (1.0 + alpha1 * rs);

             G = pre_fac * ln_res;

             double dlnd_rs = 2.0 * A * (0.5 * beta1 / rs_1 + beta2 + 1.5 * beta3 * rs_1 + (p + 1.0) * beta4 * rs_p);
             dG_drs = -2.0 * A * alpha1 * ln_res + pre_fac / ln_arg * (-dlnd_rs / (ln_denom * ln_denom));

       } // evalPW91_correlation_G 


       inline void ExchangeCorrelation::evalPW91_correlation_G1(const double& rs, double& G, double& dG_drs, const double& p){

            const double PW91_A      = 0.031091;
            const double PW91_alpha1 = 0.21370;
            const double PW91_beta1  = 7.5957;
            const double PW91_beta2  = 3.5876;
            const double PW91_beta3  = 1.6382;
            const double PW91_beta4  = 0.49294;
            
            evalPW91_correlation_G(rs,G,dG_drs,PW91_A,PW91_alpha1,PW91_beta1,PW91_beta2,PW91_beta3,PW91_beta4,p);
           
           
       }

       inline void ExchangeCorrelation::evalPW91_correlation_G2(const double& rs, double& G, double& dG_drs, const double& p){

            const double PW91_A      = 0.015545;
            const double PW91_alpha1 = 0.20548;
            const double PW91_beta1  = 14.1189;
            const double PW91_beta2  = 6.1977;
            const double PW91_beta3  = 3.3662;
            const double PW91_beta4  = 0.62517;
            
            evalPW91_correlation_G(rs,G,dG_drs,PW91_A,PW91_alpha1,PW91_beta1,PW91_beta2,PW91_beta3,PW91_beta4,p);
           
           
       }

       inline void ExchangeCorrelation::evalPW91_correlation_G3(const double& rs, double& G, double& dG_drs, const double& p){

            const double PW91_A      = 0.016887;
            const double PW91_alpha1 = 0.11125;
            const double PW91_beta1  = 10.357;
            const double PW91_beta2  = 3.6231;
            const double PW91_beta3  = 0.88026;
            const double PW91_beta4  = 0.49671;
            
            evalPW91_correlation_G(rs,G,dG_drs,PW91_A,PW91_alpha1,PW91_beta1,PW91_beta2,PW91_beta3,PW91_beta4,p);
           
           
       }



        // PW91 local-density parametrization of exchange

        inline void ExchangeCorrelation::evalPW91_exchange(const double& rho, double& fx, double& dfx_drho) {

            const double ex_par = -0.9163305865662857998 * 0.5; // in Hartree


            if (rho <= 0.0) {
                return;
            } else {
                // LDA exchange with slightly different prefactor
                evalLDA_exchange(rho, fx, dfx_drho, ex_par);

            }
        } // evalPW91_exchange





        /* ======= LOCAL DENSITY APPROXIMATION ======== */


        // LDA exchange and correlation separate contributions

        inline void ExchangeCorrelation::evalLDA(const double& rho, double& f, double& df_drho) {

            if (rho <= 0.0) {
                return;
            } else {

                // get exchange part
                double fx = 0.0;
                double dfx_drho = 0.0;
                //const double ex_par = -0.9164;
                const double ex_par = -0.4582; // in Hartree
                evalLDA_exchange(rho, fx, dfx_drho, ex_par);

                // get correlation part
                double fc = 0.0;
                double dfc_drho = 0.0;
                evalLDA_correlation(rho, fc, dfc_drho);

                f       = fx + fc;
                df_drho = dfx_drho + dfc_drho;

            }
        }




        // LDA CORRELATION

        inline void ExchangeCorrelation::evalLDA_exchange(const double& rho, double& fx, double& dfx_drho, const double& ex_par) {

            if (rho <= 0.0) {
                return;
            } else {
                double rs;
                double drs_dr;
                set_rs(rho, rs, drs_dr);

                double ex = ex_par / rs;
                double dex_drs = -ex_par / (rs * rs);
                double dex_dr = dex_drs * drs_dr;

                fx = rho * ex;
                dfx_drho = ex + rho * dex_dr;
            }

        }
        
        // LDA EXCHANGE

        inline void ExchangeCorrelation::evalLDA_correlation(const double& rho, double& fc, double& dfc_drho) {


            if (rho <= 0.0) {
                return;
            } else {

                //const double Hartree = 2.0;
                const double ec1_p = -0.1423; // * Hartree PZ_gamma_U
                const double ec1_sq = 1.0529; // PZ_beta1_U
                const double ec1_lin = 0.3334; // PX_beta2_U
                const double ec2_p = -0.0480; // * Hartree PZ_B_U
                const double ec2_log = 0.0311; // * Hartree PZ_A_U
                const double ec2_lin = -0.0116; // * Hartree PZ_D_U
                const double ec2_lig = 0.0020; // * Hartree PZ_C_U

                double ec;
                double dec_drs;

                double rs;
                double drs_dr;
                set_rs(rho, rs, drs_dr);

                if (rs >= 1.0) {

                    double sqrt_rs = sqrt(rs);
                    double ec_denom = 1.0 / (1.0 + ec1_sq * sqrt_rs + ec1_lin * rs);
                    ec = ec1_p * ec_denom;
                    dec_drs = -ec1_p * ec_denom * ec_denom * (0.5 * ec1_sq / sqrt_rs + ec1_lin);

                } else {

                    double ln_rs = log(rs);
                    ec = ec2_p + ec2_log * ln_rs + ec2_lin * rs + ec2_lig * rs * ln_rs;
                    dec_drs = ec2_log / rs + ec2_lin + ec2_lig * (1.0 + ln_rs);

                }


                double dec_drho = dec_drs * drs_dr;
                fc = rho * ec;
                dfc_drho = ec + rho * dec_drho;


            }

        }



        // END OF LOCAL DENSITY APPROXIMATION

        /*    A SET OF HELPER FUNCTIONS FOR EVALUATION OF XC */
   
        
        inline void ExchangeCorrelation::set_t(const double& rho, double& gphi, double& nrho, double& t, double& dt_drho, double& dt_dgphi, double& dt_dnrho  ){
            
            const double t_cf = 1.007715881368979474 / 4.0 ;
            const double c7_6 = 7.0/6.0;
            
            t = t_cf * nrho / pow(rho,c7_6) / gphi;
            dt_drho = -c7_6 * t / rho;
            dt_dgphi = -t / gphi;
            dt_dnrho = t / nrho;
            
            
        }
        
        
        inline void ExchangeCorrelation::set_gphi(const double& zeta, double& gphi, double& dgphi_dzeta){
  
            const double c2_3 = 2.0/3.0;
            const double c1_3 = 1.0/3.0;
            
            gphi = 0.5 * ( pow(1.0+zeta,c2_3) + pow(1.0-zeta,c2_3) );
            if ( std::abs(zeta) == 1.0 ){
                // derivative is infinite, set to large value
                if ( zeta < 0.0 ) dgphi_dzeta = -1.0e10;
                if ( zeta > 0.0 ) dgphi_dzeta = 1.0e10;
                
            } else {
                
                dgphi_dzeta = (pow(1.0+zeta,-c1_3) - pow(1.0-zeta,-c1_3) ) / 3.0;
                
            }
            
            
        } // set_gphi

        inline void ExchangeCorrelation::set_fzeta(const double& zeta, double& fzeta, double& dfzeta_dzeta) {

            const double c4_3 = 4.0 / 3.0;
            const double c1_3 = 1.0 / 3.0;


            double prefac = 1.0 / (pow(2.0, c4_3) - 2.0);
            fzeta = prefac * (pow((1.0 + zeta), c4_3) + pow(1.0 - zeta, c4_3) - 2.0);
            dfzeta_dzeta = prefac * c4_3 * (pow(1.0 + zeta, c1_3) - pow(1.0 - zeta, c1_3));

        }

        inline void ExchangeCorrelation::set_s(const double& density, double& gradNorm, double& s, double& ds_drho, double& ds_dnrho) {

            const double pisq3_cubroot = 3.093667726280135977;
            s = 0.0;


            if (gradNorm > 0) {

                s = 0.5 / pisq3_cubroot * gradNorm / pow(density, 4.0 / 3.0);
                ds_drho = -4.0 / 3.0 * s / density;
                ds_dnrho = s / gradNorm;


            } else {

                s = 0.0;
                ds_drho = 0.0;
                ds_dnrho = 0.0;

            }


        }

        inline void ExchangeCorrelation::set_rs(const double& rho, double& rs, double& drs_drho) {
            const double onethird = 1.0 / 3.0;
            const double pi = boost::math::constants::pi<double>();


            if (rho <= 0.0) {

                cerr << " Vanishing or negative density!" << endl;
                exit(1);

            } else {
                rs = pow(3.0 / (4.0 * pi * rho), onethird);
                drs_drho = -rs / rho / 3.0;
            }


        }

        void ExchangeCorrelation::setNormGradient(const double& drho_dX, const double& drho_dY, const double& drho_dZ, double& norm, double& drho_dX_normal, double& drho_dY_normal, double& drho_dZ_normal) {

            // get norm of gradient
            norm = sqrt(drho_dX * drho_dX + drho_dY * drho_dY + drho_dZ * drho_dZ);

            if (norm > 0.0) {
                // get normalized gradient
                drho_dX_normal = drho_dX / norm;
                drho_dY_normal = drho_dY / norm;
                drho_dZ_normal = drho_dZ / norm;
            }

        }


    }
}
#endif	/* EXCHANGE_CORRELATION_H */
