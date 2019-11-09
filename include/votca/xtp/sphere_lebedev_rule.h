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

#pragma once
#ifndef VOTCA_XTP_LEBEDEV_H
#define VOTCA_XTP_LEBEDEV_H

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
      const QMMolecule &atoms, const std::string &type) const;
  GridContainers::spherical_grid CalculateUnitSphereGrid(
      const std::string &element, const std::string &type) const;
  GridContainers::spherical_grid CalculateUnitSphereGrid(Index order) const;

  Index Type2MaxOrder(const std::string &element,
                      const std::string &type) const;

  Index getIndexFromOrder(Index order) const {
    if (Order2Index.count(order)) {
      return Order2Index.at(order);
    } else {
      throw std::runtime_error("No Index for Order " + std::to_string(order));
    }
  }
  Index getOrderFromIndex(Index index) const {
    if (Index2Order.count(index)) {
      return Index2Order.at(index);
    } else {
      throw std::runtime_error("No Order for Index " + std::to_string(index));
    }
  }

 private:
  Index Type2MaxOrder(const std::map<std::string, Index> &map,
                      const std::string &element) const;

  Index available_table(Index rule) const;
  Index gen_oh(Index code, double a, double b, double v, double *x, double *y,
               double *z, double *w) const;
  Eigen::Matrix4Xd ld_by_order(Index order) const;
  void ld0006(double *x, double *y, double *z, double *w) const;
  void ld0014(double *x, double *y, double *z, double *w) const;
  void ld0026(double *x, double *y, double *z, double *w) const;
  void ld0038(double *x, double *y, double *z, double *w) const;
  void ld0050(double *x, double *y, double *z, double *w) const;
  void ld0074(double *x, double *y, double *z, double *w) const;
  void ld0086(double *x, double *y, double *z, double *w) const;
  void ld0110(double *x, double *y, double *z, double *w) const;
  void ld0146(double *x, double *y, double *z, double *w) const;
  void ld0170(double *x, double *y, double *z, double *w) const;
  void ld0194(double *x, double *y, double *z, double *w) const;
  void ld0230(double *x, double *y, double *z, double *w) const;
  void ld0266(double *x, double *y, double *z, double *w) const;
  void ld0302(double *x, double *y, double *z, double *w) const;
  void ld0350(double *x, double *y, double *z, double *w) const;
  void ld0434(double *x, double *y, double *z, double *w) const;
  void ld0590(double *x, double *y, double *z, double *w) const;
  void ld0770(double *x, double *y, double *z, double *w) const;
  void ld0974(double *x, double *y, double *z, double *w) const;
  void ld1202(double *x, double *y, double *z, double *w) const;
  void ld1454(double *x, double *y, double *z, double *w) const;
  void ld1730(double *x, double *y, double *z, double *w) const;
  void ld2030(double *x, double *y, double *z, double *w) const;
  void ld2354(double *x, double *y, double *z, double *w) const;
  void ld2702(double *x, double *y, double *z, double *w) const;
  void ld3074(double *x, double *y, double *z, double *w) const;
  void ld3470(double *x, double *y, double *z, double *w) const;
  void ld3890(double *x, double *y, double *z, double *w) const;
  void ld4334(double *x, double *y, double *z, double *w) const;
  void ld4802(double *x, double *y, double *z, double *w) const;
  void ld5294(double *x, double *y, double *z, double *w) const;
  void ld5810(double *x, double *y, double *z, double *w) const;
  Index precision_table(Index rule) const;
  Index order_table(Index rule) const;
  Eigen::Vector2d Cartesian2SphericalAngle(
      const Eigen::Vector3d &r) const;  // phi=Vector[0] theta=Vector[1]

  std::map<std::string, Index> MediumOrder;
  std::map<std::string, Index> CoarseOrder;
  std::map<std::string, Index> XcoarseOrder;
  std::map<std::string, Index> FineOrder;
  std::map<std::string, Index> XfineOrder;
  std::map<Index, Index> Order2Index;
  std::map<Index, Index> Index2Order;

  void FillOrder2Index();

  void FillIndex2Order();

  void FillOrders();

  void FillMediumOrder();

  void FillFineOrder();

  void FillXfineOrder();

  void FillCoarseOrder();

  void FillXcoarseOrder();
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_LEBEDEV_H
