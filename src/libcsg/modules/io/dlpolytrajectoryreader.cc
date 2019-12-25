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

#include "dlpolytrajectoryreader.h"
#include "../../../../include/votca/csg/boundarycondition.h"
#include "../../../../include/votca/csg/topology.h"
#include <boost/filesystem/convenience.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <votca/tools/constants.h>
#include <votca/tools/getline.h>

namespace votca {
namespace csg {

using namespace std;

bool DLPOLYTrajectoryReader::Open(const string &file)
// open the original dlpoly configuration or trajectory file
// NOTE: allowed file naming - <name>.dlpc or <name>.dlph (convention:
// ".dlpc"="CONFIG", ".dlph"="HISTORY")
{
  boost::filesystem::path filepath(file.c_str());
  string inp_name = "HISTORY";

  if (boost::filesystem::extension(filepath).size() == 0) {

    throw std::ios_base::failure(
        "Error on opening dlpoly file '" + file +
        "' - extension is expected, use .dlph or .dlpc");

  } else if (boost::filesystem::extension(filepath) == ".dlpc") {

    _isConfig = true;
    inp_name = "CONFIG";

  } else if (boost::filesystem::extension(filepath) == ".dlph") {

    _isConfig = false;

  } else {
    throw std::ios_base::failure("Error on opening dlpoly file '" + file +
                                 "' - wrong extension, use .dlph or .dlpc");
  }

  if (boost::filesystem::basename(filepath).size() == 0) {
    if (filepath.parent_path().string().size() == 0) {
      _fname = inp_name;
    } else {
      _fname = filepath.parent_path().string() + "/" + inp_name;
    }
  } else {
    _fname = file;
  }

  _fl.open(_fname);
  if (!_fl.is_open()) {
    throw std::ios_base::failure("Error on opening dlpoly file '" + _fname +
                                 "'");
  }
  return true;
}

void DLPOLYTrajectoryReader::Close() { _fl.close(); }

bool DLPOLYTrajectoryReader::FirstFrame(Topology &conf) {
  _first_frame = true;
  bool res = NextFrame(conf);
  _first_frame = false;
  return res;
}

bool DLPOLYTrajectoryReader::NextFrame(Topology &conf) {
  static bool hasVs = false;
  static bool hasFs = false;
  static Index mavecs =
      0;  // number of 3d vectors per atom = keytrj in DL_POLY manuals
  static Index mpbct = 0;   // cell PBC type = imcon in DL_POLY manuals
  static Index matoms = 0;  // number of atoms/beads in a frame
  const double scale = tools::conv::ang2nm;

  static Index nerrt = 0;

  string line;

  BoundaryCondition::eBoxtype pbc_type = BoundaryCondition::typeAuto;

  if (_first_frame) {

    getline(_fl, line);  // title

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : '" << line
         << "' - header" << endl;
#endif

    getline(_fl, line);  // 2nd header line

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : '" << line
         << "' - directives line" << endl;
#endif

    tools::Tokenizer tok(line, " \t");
    vector<string> fields = tok.ToVector();

    if (fields.size() < 3) {
      throw std::runtime_error("Error: too few directive switches (<3) in '" +
                               _fname + "' header (check its 2-nd line)");
    }

    mavecs = boost::lexical_cast<Index>(fields[0]);
    mpbct = boost::lexical_cast<Index>(fields[1]);
    matoms = boost::lexical_cast<Index>(fields[2]);

    hasVs = (mavecs > 0);  // 1 or 2 => in DL_POLY frame velocity vector follows
                           // coords for each atom/bead
    hasFs = (mavecs > 1);  // 2      => in DL_POLY frame force vector follows
                           // velocities for each atom/bead

#ifdef DEBUG
    if (hasVs != conf.HasVel() || hasFs != conf.HasForce()) {
      cout << "WARNING: N of atom vectors (keytrj) in '" << _fname
           << "' header differs from that read with topology" << endl;
    }
#endif

    conf.SetHasVel(hasVs);
    conf.SetHasForce(hasFs);

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : keytrj - " << mavecs
         << ", hasV - " << conf.HasVel() << ", hasF - " << conf.HasForce()
         << endl;
#endif

    if (matoms != conf.BeadCount()) {
      throw std::runtime_error("Number of atoms/beads in '" + _fname +
                               "' header differs from that read with topology");
    }

    if (mpbct == 0) {
      pbc_type = BoundaryCondition::typeOpen;
    } else if (mpbct == 1 || mpbct == 2) {
      pbc_type = BoundaryCondition::typeOrthorhombic;
    } else if (mpbct == 3) {
      pbc_type = BoundaryCondition::typeTriclinic;
    }

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : pbc_type (imcon) - '"
         << pbc_type << "'" << endl;

    if (pbc_type != conf.getBoxType())
      cout << "WARNING: PBC type in dlpoly file '" << _fname
           << "' header differs from that read with topology" << endl;
// throw std::runtime_error("Error: Boundary conditions in '"+_fname+"'
// header differs from that read with topology");
#endif
  } else if (_isConfig) {

    return false;
  }
  // read normal frame

  if (!_isConfig) {
    getline(_fl, line);  // timestep line - only present in HISTORY, and not in
                         // CONFIG
#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : '" << line << "'"
         << endl;
#endif
  }

  if (!_fl.eof()) {
    Index nstep;
    Index natoms;
    Index navecs;

    if (_isConfig) {
      // use the above read specs from the header, and skip the data missing in
      // CONFIG

      natoms = matoms;
      navecs = mavecs;

      conf.SetHasVel(hasVs);
      conf.SetHasForce(hasFs);

#ifdef DEBUG
      cout << "Read from CONFIG: traj_key - " << navecs << ", hasV - "
           << conf.HasVel() << ", hasF - " << conf.HasForce() << endl;
#endif

    } else {

      tools::Tokenizer tok(line, " \t");
      vector<string> fields = tok.ToVector();

      if (fields.size() < 6) {
        throw std::runtime_error(
            "Error: too few directive switches (<6) in 'timestep' record");
      }

      nstep = boost::lexical_cast<Index>(fields[1]);
      natoms = boost::lexical_cast<Index>(fields[2]);
      navecs = boost::lexical_cast<Index>(fields[3]);
      Index npbct = boost::lexical_cast<Index>(fields[4]);
      double dtime =
          stod(fields[5]);  // normally it is the 5-th column in 'timestep' line
      double stime =
          stod(fields[fields.size() - 1]);  // normally it is the last
                                            // column in 'timestep' line

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : natoms = " << natoms
           << ", levcfg = " << fields[3];
      cout << ", dt = " << fields[5] << ", time = " << stime << endl;
#endif

      if (natoms != conf.BeadCount()) {
        throw std::runtime_error(
            "Error: N of atoms/beads in '" + _fname +
            "' header differs from that found in topology");
      }
      if (natoms != matoms) {
        throw std::runtime_error(
            "Error: N of atoms/beads in '" + _fname +
            "' header differs from that found in the frame");
      }
      if (navecs != mavecs) {
        throw std::runtime_error(
            "Error: N of atom vectors (keytrj) in '" + _fname +
            "' header differs from that found in the frame");
      }
      if (npbct != mpbct) {
        throw std::runtime_error(
            "Error: boundary conditions (imcon) in '" + _fname +
            "' header differs from that found in the frame");
      }

      // total time - calculated as product due to differences between DL_POLY
      // versions in HISTORY formats
      conf.setTime(double(nstep) * dtime);
      conf.setStep(nstep);

      if (std::abs(stime - conf.getTime()) > 1.e-8) {
        nerrt++;
        if (nerrt < 11) {
          cout << "Check: nstep = " << nstep << ", dt = " << dtime
               << ", time = " << stime << " (correct?)" << endl;
        } else if (nerrt == 11) {
          cout << "Check: timestep - more than 10 mismatches in total time "
                  "found..."
               << endl;
        }
      }

      if (npbct == 0) {
        pbc_type = BoundaryCondition::typeOpen;
      } else if (npbct == 1 || npbct == 2) {
        pbc_type = BoundaryCondition::typeOrthorhombic;
      } else if (npbct == 3) {
        pbc_type = BoundaryCondition::typeTriclinic;
      }
    }

    Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
    for (Index i = 0; i < 3; i++) {  // read 3 box/cell lines

      getline(_fl, line);

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << line
           << "' - box vector # " << i + 1 << endl;
#endif

      if (_fl.eof()) {
        throw std::runtime_error("Error: unexpected EOF in dlpoly file '" +
                                 _fname + "', when reading box vector" +
                                 boost::lexical_cast<string>(i));
      }

      tools::Tokenizer tok(line, " \t");
      vector<double> fields;
      tok.ConvertToVector<double>(fields);
      // Angs -> nm
      box.col(i) = scale * Eigen::Vector3d(fields[0], fields[1], fields[2]);
    }

    conf.setBox(box, pbc_type);

    for (Index i = 0; i < natoms; i++) {

      {
        getline(_fl, line);  // atom header line

#ifdef DEBUG
        cout << "Read from dlpoly file '" << _fname << "' : '" << line << "'"
             << endl;
#endif

        if (_fl.eof()) {
          throw std::runtime_error("Error: unexpected EOF in dlpoly file '" +
                                   _fname + "', when reading atom/bead # " +
                                   boost::lexical_cast<string>(i + 1));
        }

        tools::Tokenizer tok(line, " \t");
        vector<string> fields = tok.ToVector();
        Index id = boost::lexical_cast<Index>(fields[1]);
        if (i + 1 != id) {
          throw std::runtime_error(
              "Error: unexpected atom/bead index in dlpoly file '" + _fname +
              "' : expected " + boost::lexical_cast<string>(i + 1) +
              " but got " + boost::lexical_cast<string>(id));
        }
      }

      Bead *b = conf.getBead(i);
      Eigen::Matrix3d atom_vecs = Eigen::Matrix3d::Zero();
      for (Index j = 0; j < std::min(navecs, Index(2)) + 1; j++) {

        getline(_fl, line);  // read atom positions

#ifdef DEBUG
        cout << "Read from dlpoly file '" << _fname << "' : '" << line << "'"
             << endl;
#endif

        if (_fl.eof()) {
          throw std::runtime_error(
              "Error: unexpected EOF in dlpoly file '" + _fname +
              "', when reading atom/bead vector # " +
              boost::lexical_cast<string>(j) + " of atom " +
              boost::lexical_cast<string>(i + 1));
        }

        vector<double> fields;
        tools::Tokenizer tok(line, " \t");
        tok.ConvertToVector<double>(fields);
        // Angs -> nm
        atom_vecs.col(j) =
            scale * Eigen::Vector3d(fields[0], fields[1], fields[2]);
      }

      b->setPos(atom_vecs.col(0));
#ifdef DEBUG
      cout << "Crds from dlpoly file '" << _fname << "' : " << atom_vecs.col(0)
           << endl;
#endif
      if (navecs > 0) {
        b->setVel(atom_vecs.col(1));
#ifdef DEBUG
        cout << "Vels from dlpoly file '" << _fname
             << "' : " << atom_vecs.col(1) << endl;
#endif
        if (navecs > 1) {
          b->setF(atom_vecs.col(2));
#ifdef DEBUG
          cout << "Frcs from dlpoly file '" << _fname
               << "' : " << atom_vecs.col(2) << endl;
#endif
        }
      }
    }
  }

  return !_fl.eof();
}

}  // namespace csg
}  // namespace votca
