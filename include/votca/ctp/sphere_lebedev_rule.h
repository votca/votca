/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __CTP_LEBEDEV__H
#define	__CTP_LEBEDEV__H


#include <votca/tools/property.h>
#include <votca/ctp/grid_containers.h>
#include <boost/math/constants/constants.hpp>
using namespace std;


namespace votca {
    namespace ctp {

        class LebedevGrid {
        public:

            LebedevGrid() {
                FillOrders();
            };

            void getSphericalGrid(vector<QMAtom* > _atoms, string type, GridContainers& _grids);
            void getUnitSphereGrid(string element, string type, std::vector<double>& _theta, std::vector<double>& _phi, std::vector<double>& _weight);



        private:
            const double pi = boost::math::constants::pi<double>();
            int available_table(int rule);
            int gen_oh(int code, double a, double b, double v, double *x,
                    double *y, double *z, double *w);
            void ld_by_order(int order, double *x, double *y, double *z, double *w);
            void ld0006(double *x, double *y, double *z, double *w);
            void ld0014(double *x, double *y, double *z, double *w);
            void ld0026(double *x, double *y, double *z, double *w);
            void ld0038(double *x, double *y, double *z, double *w);
            void ld0050(double *x, double *y, double *z, double *w);
            void ld0074(double *x, double *y, double *z, double *w);
            void ld0086(double *x, double *y, double *z, double *w);
            void ld0110(double *x, double *y, double *z, double *w);
            void ld0146(double *x, double *y, double *z, double *w);
            void ld0170(double *x, double *y, double *z, double *w);
            void ld0194(double *x, double *y, double *z, double *w);
            void ld0230(double *x, double *y, double *z, double *w);
            void ld0266(double *x, double *y, double *z, double *w);
            void ld0302(double *x, double *y, double *z, double *w);
            void ld0350(double *x, double *y, double *z, double *w);
            void ld0434(double *x, double *y, double *z, double *w);
            void ld0590(double *x, double *y, double *z, double *w);
            void ld0770(double *x, double *y, double *z, double *w);
            void ld0974(double *x, double *y, double *z, double *w);
            void ld1202(double *x, double *y, double *z, double *w);
            void ld1454(double *x, double *y, double *z, double *w);
            void ld1730(double *x, double *y, double *z, double *w);
            void ld2030(double *x, double *y, double *z, double *w);
            void ld2354(double *x, double *y, double *z, double *w);
            void ld2702(double *x, double *y, double *z, double *w);
            void ld3074(double *x, double *y, double *z, double *w);
            void ld3470(double *x, double *y, double *z, double *w);
            void ld3890(double *x, double *y, double *z, double *w);
            void ld4334(double *x, double *y, double *z, double *w);
            void ld4802(double *x, double *y, double *z, double *w);
            void ld5294(double *x, double *y, double *z, double *w);
            void ld5810(double *x, double *y, double *z, double *w);
            int order_table(int rule);
            int precision_table(int rule);
            void timestamp();
            void xyz_to_tp(double x, double y, double z, double *t, double *p);

            int getOrder(string element, string type);



            std::map<std::string, int> MediumOrder;

            inline void FillOrders() {

                FillMediumOrder();

            }

            inline void FillMediumOrder() {

                // order for H, He (not given in NWChem, assuming same as 1st row)
                MediumOrder["H"] = 434;
                MediumOrder["He"] = 434;

                // orders for 1st row elements taken from NWChem
                MediumOrder["Li"] = 434;
                MediumOrder["Be"] = 434;
                MediumOrder["B"] = 434;
                MediumOrder["C"] = 434;
                MediumOrder["N"] = 434;
                MediumOrder["O"] = 434;
                MediumOrder["F"] = 434;
                MediumOrder["Ne"] = 434;

                // orders for 2nd row elements taken from NWChem
                MediumOrder["Na"] = 434;
                MediumOrder["Mg"] = 434;
                MediumOrder["Al"] = 434;
                MediumOrder["Si"] = 434;
                MediumOrder["P"] = 434;
                MediumOrder["S"] = 434;
                MediumOrder["Cl"] = 434;
                MediumOrder["Ar"] = 434;


                // orders for 3rd row elements taken from NWChem
                MediumOrder["K"] = 590;
                MediumOrder["Ca"] = 590;
                MediumOrder["Sc"] = 590;
                MediumOrder["Ti"] = 590;
                MediumOrder["V"] = 590;
                MediumOrder["Cr"] = 590;
                MediumOrder["Mn"] = 590;
                MediumOrder["Fe"] = 590;
                MediumOrder["Co"] = 590;
                MediumOrder["Ni"] = 590;
                MediumOrder["Cu"] = 590;
                MediumOrder["Zn"] = 590;
                MediumOrder["Ga"] = 590;
                MediumOrder["Ge"] = 590;
                MediumOrder["As"] = 590;
                MediumOrder["Se"] = 590;
                MediumOrder["Br"] = 590;
                MediumOrder["Kr"] = 590;


            }



        };

    }
}
#endif	/* LEBEDEV_H */