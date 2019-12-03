/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _TRAJ_FORCE_H
#define _TRAJ_FORCE_H

#include <boost/numeric/ublas/vector.hpp>
#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/property.h>

using namespace votca::csg;
using namespace std;

/**
   \brief Adds/subtracts forces from given atomistic trajectories
**/

class TrajForce : public CsgApplication {
 public:
  string ProgramName() override { return "traj_force"; }
  void HelpText(ostream &out) override {
    out << "Adds/subtracts forces from given atomistic trajectories";
  }

  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return false; }

  void Initialize(void) override;
  bool EvaluateOptions() override;

  /// \brief called before the first frame
  void BeginEvaluate(Topology *top, Topology *top_atom) override;
  /// \brief called after the last frame
  void EndEvaluate() override;
  /// \brief called for each frame which is mapped
  void EvalConfiguration(Topology *conf, Topology *conf_atom) override;

 protected:
  /// \brief Scaling of forces, +1 for addition and -1 for subtraction
  double _scale;
  /// \brief Write results to output files
  void WriteOutFiles();

  void OpenForcesTrajectory();
  Topology _top_force;
  TrajectoryReader *_trjreader_force;
  TrajectoryWriter *_trjwriter;
};

#endif /* _TRAJ_FORCE_H */
