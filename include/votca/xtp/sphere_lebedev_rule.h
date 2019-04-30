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

#ifndef __XTP_LEBEDEV__H
#define __XTP_LEBEDEV__H

#include <boost/math/constants/constants.hpp>
#include <votca/tools/property.h>
#include <votca/xtp/grid_containers.h>
#include <votca/xtp/qmatom.h>

namespace votca {
namespace xtp {

class LebedevGrid {
 public:
  LebedevGrid() {
    FillOrders();
    FillOrder2Index();
    FillIndex2Order();
  };

  std::map<std::string, GridContainers::spherical_grid> CalculateSphericalGrids(
      std::vector<QMAtom *> atoms, const std::string &type);
  GridContainers::spherical_grid CalculateUnitSphereGrid(
      const std::string &element, const std::string &type);
  GridContainers::spherical_grid CalculateUnitSphereGrid(int order);

  int Type2MaxOrder(const std::string &element, const std::string &type);
  int getIndexFromOrder(int order) { return Order2Index.at(order); }
  int getOrderFromIndex(int index) { return Index2Order.at(index); }

 private:
  int available_table(int rule);
  int gen_oh(int code, double a, double b, double v, double *x, double *y,
             double *z, double *w);
  Eigen::Matrix4Xd ld_by_order(int order);
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
  int precision_table(int rule);
  int order_table(int rule);
  Eigen::Vector2d Cartesian2SphericalAngle(
      const Eigen::Vector3d &r);  // phi=Vector[0] theta=Vector[1]

  std::map<std::string, int> MediumOrder;
  std::map<std::string, int> CoarseOrder;
  std::map<std::string, int> XcoarseOrder;
  std::map<std::string, int> FineOrder;
  std::map<std::string, int> XfineOrder;
  std::map<int, int> Order2Index;
  std::map<int, int> Index2Order;

  inline void FillOrder2Index() {
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

  inline void FillIndex2Order() {
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
    FillCoarseOrder();
    FillXcoarseOrder();
    FillFineOrder();
    FillXfineOrder();
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

    // 4th row (selection)
    MediumOrder["Ag"] = 590;
  }

  inline void FillFineOrder() {
    // order for H, He (not given in NWChem, assuming same as 1st row)
    FineOrder["H"] = 590;
    FineOrder["He"] = 590;

    // orders for 1st row elements taken from NWChem
    FineOrder["Li"] = 590;
    FineOrder["Be"] = 590;
    FineOrder["B"] = 590;
    FineOrder["C"] = 590;
    FineOrder["N"] = 590;
    FineOrder["O"] = 590;
    FineOrder["F"] = 590;
    FineOrder["Ne"] = 590;

    // orders for 2nd row elements taken from NWChem
    FineOrder["Na"] = 770;
    FineOrder["Mg"] = 770;
    FineOrder["Al"] = 770;
    FineOrder["Si"] = 770;
    FineOrder["P"] = 770;
    FineOrder["S"] = 770;
    FineOrder["Cl"] = 770;
    FineOrder["Ar"] = 770;

    // orders for 3rd row elements taken from NWChem
    FineOrder["K"] = 974;
    FineOrder["Ca"] = 974;
    FineOrder["Sc"] = 974;
    FineOrder["Ti"] = 974;
    FineOrder["V"] = 974;
    FineOrder["Cr"] = 974;
    FineOrder["Mn"] = 974;
    FineOrder["Fe"] = 974;
    FineOrder["Co"] = 974;
    FineOrder["Ni"] = 974;
    FineOrder["Cu"] = 974;
    FineOrder["Zn"] = 974;
    FineOrder["Ga"] = 974;
    FineOrder["Ge"] = 974;
    FineOrder["As"] = 974;
    FineOrder["Se"] = 974;
    FineOrder["Br"] = 974;
    FineOrder["Kr"] = 974;

    // 4th row
    FineOrder["Ag"] = 974;
  }
  inline void FillXfineOrder() {
    // order for H, He (not given in NWChem, assuming same as 1st row)
    XfineOrder["H"] = 1202;
    XfineOrder["He"] = 1202;

    // orders for 1st row elements taken from NWChem
    XfineOrder["Li"] = 1202;
    XfineOrder["Be"] = 1202;
    XfineOrder["B"] = 1202;
    XfineOrder["C"] = 1202;
    XfineOrder["N"] = 1202;
    XfineOrder["O"] = 1202;
    XfineOrder["F"] = 1202;
    XfineOrder["Ne"] = 1202;

    // orders for 2nd row elements taken from NWChem
    XfineOrder["Na"] = 1454;
    XfineOrder["Mg"] = 1454;
    XfineOrder["Al"] = 1454;
    XfineOrder["Si"] = 1454;
    XfineOrder["P"] = 1454;
    XfineOrder["S"] = 1454;
    XfineOrder["Cl"] = 1454;
    XfineOrder["Ar"] = 1454;

    // orders for 3rd row elements taken from NWChem
    XfineOrder["K"] = 1454;
    XfineOrder["Ca"] = 1454;
    XfineOrder["Sc"] = 1454;
    XfineOrder["Ti"] = 1454;
    XfineOrder["V"] = 1454;
    XfineOrder["Cr"] = 1454;
    XfineOrder["Mn"] = 1454;
    XfineOrder["Fe"] = 1454;
    XfineOrder["Co"] = 1454;
    XfineOrder["Ni"] = 1454;
    XfineOrder["Cu"] = 1454;
    XfineOrder["Zn"] = 1454;
    XfineOrder["Ga"] = 1454;
    XfineOrder["Ge"] = 1454;
    XfineOrder["As"] = 1454;
    XfineOrder["Se"] = 1454;
    XfineOrder["Br"] = 1454;
    XfineOrder["Kr"] = 1454;

    // 4th row
    XfineOrder["Ag"] = 1454;
  }

  inline void FillCoarseOrder() {
    // order for H, He (not given in NWChem, assuming same as 1st row)
    CoarseOrder["H"] = 302;
    CoarseOrder["He"] = 302;

    // orders for 1st row elements taken from NWChem
    CoarseOrder["Li"] = 302;
    CoarseOrder["Be"] = 302;
    CoarseOrder["B"] = 302;
    CoarseOrder["C"] = 302;
    CoarseOrder["N"] = 302;
    CoarseOrder["O"] = 302;
    CoarseOrder["F"] = 302;
    CoarseOrder["Ne"] = 302;

    // orders for 2nd row elements taken from NWChem
    CoarseOrder["Na"] = 302;
    CoarseOrder["Mg"] = 302;
    CoarseOrder["Al"] = 302;
    CoarseOrder["Si"] = 302;
    CoarseOrder["P"] = 302;
    CoarseOrder["S"] = 302;
    CoarseOrder["Cl"] = 302;
    CoarseOrder["Ar"] = 302;

    // orders for 3rd row elements taken from NWChem
    CoarseOrder["K"] = 302;
    CoarseOrder["Ca"] = 302;
    CoarseOrder["Sc"] = 302;
    CoarseOrder["Ti"] = 302;
    CoarseOrder["V"] = 302;
    CoarseOrder["Cr"] = 302;
    CoarseOrder["Mn"] = 302;
    CoarseOrder["Fe"] = 302;
    CoarseOrder["Co"] = 302;
    CoarseOrder["Ni"] = 302;
    CoarseOrder["Cu"] = 302;
    CoarseOrder["Zn"] = 302;
    CoarseOrder["Ga"] = 302;
    CoarseOrder["Ge"] = 302;
    CoarseOrder["As"] = 302;
    CoarseOrder["Se"] = 302;
    CoarseOrder["Br"] = 302;
    CoarseOrder["Kr"] = 302;

    // 4th row
    CoarseOrder["Ag"] = 302;
  }

  inline void FillXcoarseOrder() {
    // order for H, He (not given in NWChem, assuming same as 1st row)
    XcoarseOrder["H"] = 194;
    XcoarseOrder["He"] = 194;

    // orders for 1st row elements taken from NWChem
    XcoarseOrder["Li"] = 194;
    XcoarseOrder["Be"] = 194;
    XcoarseOrder["B"] = 194;
    XcoarseOrder["C"] = 194;
    XcoarseOrder["N"] = 194;
    XcoarseOrder["O"] = 194;
    XcoarseOrder["F"] = 194;
    XcoarseOrder["Ne"] = 194;

    // orders for 2nd row elements taken from NWChem
    XcoarseOrder["Na"] = 194;
    XcoarseOrder["Mg"] = 194;
    XcoarseOrder["Al"] = 194;
    XcoarseOrder["Si"] = 194;
    XcoarseOrder["P"] = 194;
    XcoarseOrder["S"] = 194;
    XcoarseOrder["Cl"] = 194;
    XcoarseOrder["Ar"] = 194;

    // orders for 3rd row elements taken from NWChem
    XcoarseOrder["K"] = 194;
    XcoarseOrder["Ca"] = 194;
    XcoarseOrder["Sc"] = 194;
    XcoarseOrder["Ti"] = 194;
    XcoarseOrder["V"] = 194;
    XcoarseOrder["Cr"] = 194;
    XcoarseOrder["Mn"] = 194;
    XcoarseOrder["Fe"] = 194;
    XcoarseOrder["Co"] = 194;
    XcoarseOrder["Ni"] = 194;
    XcoarseOrder["Cu"] = 194;
    XcoarseOrder["Zn"] = 194;
    XcoarseOrder["Ga"] = 194;
    XcoarseOrder["Ge"] = 194;
    XcoarseOrder["As"] = 194;
    XcoarseOrder["Se"] = 194;
    XcoarseOrder["Br"] = 194;
    XcoarseOrder["Kr"] = 194;

    // 4th row
    XcoarseOrder["Ag"] = 194;
  }
};

}  // namespace xtp
}  // namespace votca
#endif /* LEBEDEV_H */