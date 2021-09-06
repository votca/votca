/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Local VOTCA includes
#include "votca/xtp/segmentmapper.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "bgpol.h"
#include "ewald/ewd_segment.h"
#include "ewald/unitcell.h"

namespace votca {
namespace xtp {

void BGPol::ParseOptions(const tools::Property& options) {
  _mapfile = options.get(".mapfile").as<std::string>();
  ewd_options.alpha =
      (1.0 / tools::conv::nm2bohr) * options.get(".alpha").as<double>();
  ewd_options.k_cutoff =
      (1.0 / tools::conv::nm2bohr) * options.get(".k_cutoff").as<double>();
  ewd_options.r_cutoff =
      tools::conv::nm2bohr * options.get(".r_cutoff").as<double>();
  ewd_options.sharpness = options.get(".thole_sharpness").as<double>();
  std::string shape_str = options.get(".shape").as<std::string>();
  if (shape_str == "cube") {
    ewd_options.shape = Shape::cube;
  } else if (shape_str == "sphere") {
    ewd_options.shape = Shape::sphere;
  } else if (shape_str == "xyslab") {
    ewd_options.shape = Shape::xyslab;
  } else {
    throw std::runtime_error("Unknown shape in option file\n");
  }
}

bool BGPol::Evaluate(Topology& top) {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  // Make XTP_LOG behave like std::cout
  _log.setCommonPreface("... ...");
  XTP_LOG(Log::error, _log) << std::endl;

  // Map multipole and polarization data to segments
  PolarMapper polmap(_log);
  polmap.LoadMappingFile(_mapfile);
  for (const Segment& seg : top.Segments()) {
    PolarSegment mol = polmap.map(seg, SegId(seg.getId(), "n"));
    _polar_background.push_back(mol);
  }

  // Convert data to a cartesian representation
  std::vector<EwdSegment> _ewald_background;
  for (const PolarSegment& pseg : _polar_background) {
    EwdSegment eseg(pseg);
    _ewald_background.push_back(eseg);
  }

  uc_matrix_ = top.getBox();
  UnitCell unit_cell(uc_matrix_);

  // Place atoms in the simulation box
  for (EwdSegment& seg : _ewald_background) {
    for (EwdSite& site : seg) {
      site.updatePos(unit_cell.placeCoordInBox(site.getPos()));
    }
    // Should be deleted when merged and after compare with ctp is done
    seg.calcPos();
  }

  // Polarize the neutral background
  BackgroundPolarizer BgPol(_log, unit_cell, ewd_options);
  BgPol.Polarize(_ewald_background);

  // Write the result to an hdf5 file
  std::string output_file_name = "background_polarization.hdf5";
  WriteToHdf5(output_file_name, _ewald_background);

  return true;
}

void BGPol::WriteToHdf5(std::string filename,
                        std::vector<EwdSegment> background) const {
  // Write the polarization state including all segments
  CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
  CheckpointWriter a = cpf.getWriter();
  CheckpointWriter ww = a.openChild("polar_background");
  for (EwdSegment& seg : background) {
    CheckpointWriter www =
        ww.openChild("background_" + std::to_string(seg.getId()));
    seg.WriteToCpt(www);
  }

  // Write the ewald options
  CheckpointWriter w2 = a.openChild("ewald_options");
  w2(ewd_options.alpha, "alpha");
  w2(ewd_options.k_cutoff, "k_cutoff");
  w2(ewd_options.r_cutoff, "r_cutoff");
  w2(ewd_options.sharpness, "sharpness");
  switch (ewd_options.shape) {
    case Shape::cube:
      w2(std::string("cube"), "shape");
      break;
    case Shape::xyslab:
      w2(std::string("xyslab"), "shape");
      break;
    case Shape::sphere:
      w2(std::string("sphere"), "shape");
      break;
    default:
      throw std::runtime_error("Shape not recognized!");
  }
  // Write the unit cell info
  CheckpointWriter w3 = a.openChild("unit_cell");
  w3(uc_matrix_, "unit_cell_matrix");
}

}  // namespace xtp

}  // namespace votca
