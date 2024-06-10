/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#ifndef VOTCA_XTP_REGION_H
#define VOTCA_XTP_REGION_H

// Standard includes
#include <iostream>

// VOTCA includes
#include <votca/csg/pdbwriter.h>

// Local VOTCA includes
#include "checkpoint.h"
#include "logger.h"

/**
 * \brief base class to derive regions from
 *
 *
 *
 */

namespace votca {
namespace xtp {

class QMRegion;
class PolarRegion;
class StaticRegion;

class Region {

 public:
  Region(Index id, Logger& log) : id_(id), log_(log) {};
  virtual ~Region() = default;

  virtual void WriteToCpt(CheckpointWriter& w) const = 0;

  virtual void ReadFromCpt(CheckpointReader& r) = 0;

  virtual void Initialize(const tools::Property& prop) = 0;

  virtual bool Converged() const = 0;

  virtual void Evaluate(std::vector<std::unique_ptr<Region> >& regions) = 0;

  virtual Index size() const = 0;

  virtual std::string identify() const = 0;

  virtual void WritePDB(csg::PDBWriter& writer) const = 0;

  virtual void Reset() = 0;

  virtual double charge() const = 0;

  bool Successful() const { return info_; }

  std::string ErrorMsg() const { return errormsg_; }

  void AddResults(tools::Property& prop) const;

  Index getId() const { return id_; }

  virtual double Etotal() const = 0;

  friend std::ostream& operator<<(std::ostream& out, const Region& region) {
    out << "Id: " << region.getId() << " type: " << region.identify()
        << " size: " << region.size() << " charge[e]= " << region.charge();
    return out;
  }

 protected:
  bool info_ = true;
  std::string errormsg_ = "";
  std::vector<double> ApplyInfluenceOfOtherRegions(
      std::vector<std::unique_ptr<Region> >& regions);
  virtual void AppendResult(tools::Property& prop) const = 0;
  virtual double InteractwithQMRegion(const QMRegion& region) = 0;
  virtual double InteractwithPolarRegion(const PolarRegion& region) = 0;
  virtual double InteractwithStaticRegion(const StaticRegion& region) = 0;

  Index id_ = -1;
  Logger& log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_REGION_H
