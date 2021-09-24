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
#include "votca/xtp/ewd_segment.h"
#include "votca/xtp/unitcell.h"

namespace votca {
namespace xtp {

void BGPol::ParseOptions(const tools::Property& options) {
  mapfile_ = options.get(".mapfile").as<std::string>();
  ewd_options_.alpha =
      (1.0 / tools::conv::nm2bohr) * options.get(".alpha").as<double>();
  ewd_options_.k_cutoff =
      (1.0 / tools::conv::nm2bohr) * options.get(".k_cutoff").as<double>();
  ewd_options_.r_cutoff =
      tools::conv::nm2bohr * options.get(".r_cutoff").as<double>();
  ewd_options_.sharpness = options.get(".thole_sharpness").as<double>();
  std::string shape_str = options.get(".shape").as<std::string>();
  if (shape_str == "cube") {
    ewd_options_.shape = Shape::cube;
  } else if (shape_str == "sphere") {
    ewd_options_.shape = Shape::sphere;
  } else if (shape_str == "xyslab") {
    ewd_options_.shape = Shape::xyslab;
  } else {
    throw std::runtime_error("Unknown shape in option file\n");
  }
  output_file_name_ = options.get("output_file").as<std::string>() + ".hdf5";
}

bool BGPol::Evaluate(Topology& top) {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  // Make XTPlog_ behave like std::cout
  log_.setCommonPreface("... ...");
  XTP_LOG(Log::error, log_) << std::endl;

  // Map multipole and polarization data to segments
  std::vector<PolarSegment> polar_background;
  PolarMapper polmap(log_);
  polmap.LoadMappingFile(mapfile_);
  for (const Segment& seg : top.Segments()) {
    PolarSegment mol = polmap.map(seg, SegId(seg.getId(), "n"));
    polar_background.push_back(mol);
  }

  // Polarize the neutral background
  uc_matrix_ = top.getBox();
  Background Bg(log_, top.getBox(), ewd_options_, polar_background);
  Bg.Polarize();
  Bg.writeToStateFile(output_file_name_);
  return true;
}

}  // namespace xtp

}  // namespace votca
