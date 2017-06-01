/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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



namespace votca {
    namespace xtp {

        class ExchangeCorrelation {
        public:

            ExchangeCorrelation() {
            };

            void getXC(int type, const double rho, const double drho_dX, const double drho_dY, const double drho_dZ, double& f, double& df_drho, double& df_dsigma);



        private:

            void set_rs(const double rho, double& rs, double& drs_dr);
            void set_fzeta(const double zeta, double& fzeta, double& dfzeta_dzeta);
            void set_s(const double density, double& gradNorm, double& s, double& ds_drho, double& ds_dnrho);
            void setNormGradient(const double drho_dX, const double drho_dY, const double drho_dZ, double& norm, double& drho_dX_normal, double& drho_dY_normal, double& drho_dZ_normal);
            void set_gphi(const double zeta, double& gphi, double& dgphi_dzeta);
            void set_t(const double rho, double& gphi, double& nrho, double& t, double& dt_drho, double& dt_dgphi, double& dt_dnrho  );

            void evalLDA(const double rho, double& f, double& df_drho);
            void evalLDA_exchange(const double rho, double& fx, double& dfx_drho, const double ex_par);
            void evalLDA_correlation(const double rho, double& fc, double& dfc_drho);

            void evalPW91_exchange(const double rho, double& fx, double& dfx_drho);
            void evalPW91_correlation_G1( const double rs, double& G, double& dG_drs, const double p);
            void evalPW91_correlation_G2( const double rs, double& G, double& dG_drs, const double p);
            void evalPW91_correlation_G3( const double rs, double& G, double& dG_drs, const double p);
            void evalPW91_correlation_G( const double rs, double& G, double& dG_drs, const double A, const double alpha1, const double beta1, const double beta2, const double beta3, const double beta4, const double p);
            void evalPW91_correlation_int( const double rho, const double zeta, double& ec, double& dec_drho, double& dec_dzeta );

            //void evalPBE(const double rho, const double drho_dX, const double drho_dY, const double drho_dZ, double& f, double& df_drho, double& df_drho_dX, double& df_drho_dY, double& df_drho_dZ);
            void evalPBE(const double rho, const double drho_dX, const double drho_dY, const double drho_dZ, double& f, double& df_drho, double& df_dsigma);
            void evalPBE_exchange(const double rho, const double drho_dX, const double drho_dY, const double drho_dZ, double& fx, double& dfx_drho, double& dfx_drho_dX, double& dfx_drho_dY, double& dfx_drho_dZ);
            void evalPBE_correlation(const double rho, const double drho_dX, const double drho_dY, const double drho_dZ, double& fc, double& dfc_drho, double& dfc_drho_dX,double& dfc_drho_dY,double& dfc_drho_dZ);
            void evalPBE_correlation_H(const double rho,  double& zeta,  double& nrho,  double& eps,  double& deps_drho,  double& deps_dzeta, double& H, double& dH_drho, double& dH_dzeta, double& dH_dnrho);

        };
    }
}
#endif	/* EXCHANGE_CORRELATION_H */
