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
                FillOrder2Index();
                FillIndex2Order();
            };

            void getSphericalGrid(vector<QMAtom* > _atoms, string type, GridContainers& _grids);
            void getUnitSphereGrid(string element, string type, std::vector<double>& _theta, std::vector<double>& _phi, std::vector<double>& _weight);
            void getUnitSphereGrid(int order, std::vector<double>& _theta, std::vector<double>& _phi, std::vector<double>& _weight);

            int Type2MaxOrder( string element, string type );
            int getIndexFromOrder( int order ) { return Order2Index.at(order); }
            int getOrderFromIndex( int index ) { return Index2Order.at(index); }
            

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
            std::map<int,int>          Order2Index;
            std::map<int,int>          Index2Order;

            
            inline void FillOrder2Index(){
                
                Order2Index[38] = 1;
                Order2Index[50] = 2;
                Order2Index[74] = 3;
                Order2Index[86] = 4;
                Order2Index[110] = 5;
                Order2Index[146] = 6;
                Order2Index[170] = 7;
                Order2Index[194] = 8;
                Order2Index[230] = 9;
                Order2Index[266] = 10;
                Order2Index[302] = 11;
                Order2Index[350] = 12;
                Order2Index[434] = 13;
                Order2Index[590] = 14;
                Order2Index[770] = 15;
                Order2Index[974] = 16;
                Order2Index[1202] = 17;
                Order2Index[1454] = 18;
                Order2Index[1730] = 19;
                Order2Index[2030] = 20;
                Order2Index[2354] = 21;
                Order2Index[2702] = 22;
                Order2Index[3074] = 23;
                Order2Index[3470] = 24;
                Order2Index[3890] = 25;
                Order2Index[4334] = 26;
                Order2Index[4802] = 27;
                Order2Index[5294] = 28;
                Order2Index[5810] = 29;
                
            }
            
            inline void FillIndex2Order(){
                
                Index2Order[1] = 38;
                Index2Order[2] = 50;
                Index2Order[3] = 74;
                Index2Order[4] = 86;
                Index2Order[5] = 110;
                Index2Order[6] = 146;
                Index2Order[7] = 170;
                Index2Order[8] = 194;
                Index2Order[9] = 230;
                Index2Order[10] = 266;
                Index2Order[11] = 302;
                Index2Order[12] = 350;
                Index2Order[13] = 434;
                Index2Order[14] = 590;
                Index2Order[15] = 770;
                Index2Order[16] = 974;
                Index2Order[17] = 1202;
                Index2Order[18] = 1454;
                Index2Order[19] = 1730;
                Index2Order[20] = 2030;
                Index2Order[21] = 2354;
                Index2Order[22] = 2702;
                Index2Order[23] = 3074;
                Index2Order[24] = 3470;
                Index2Order[25] = 3890;
                Index2Order[26] = 4334;
                Index2Order[27] = 4802;
                Index2Order[28] = 5294;
                Index2Order[29] = 5810;
                
            }
            
            
            
            
            
            
            
            
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