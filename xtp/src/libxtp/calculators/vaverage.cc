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

// Local VOTCA includes
#include "votca/xtp/topology.h"
#include <votca/tools/getline.h>

// Local private VOTCA includes
#include "vaverage.h"

namespace votca {
namespace xtp {

void VAverage::ParseOptions(const tools::Property& options) {

  ratefile_ = options.get(".ratefile").as<std::string>();
  occfile_ = options.get(".occfile").as<std::string>();
  outputfile_ = options.get(".output").as<std::string>();
}

std::vector<double> VAverage::ReadOccfile(std::string filename) const {
  std::vector<double> result;
  std::ifstream intt;
  intt.open(filename);
  if (!intt.is_open()) {
    throw std::runtime_error("File:" + filename + " could not be opened");
  }
  std::string line;
  Index id = 0;
  while (intt.good()) {

    tools::getline(intt, line);
    std::vector<std::string> split = tools::Tokenizer(line, " \t").ToVector();
    if (!split.size() || split[0] == "#" || split[0].substr(0, 1) == "#") {
      continue;
    }
    if (split.size() != 2) {
      throw std::runtime_error("Row should only contain id and occupation");
    }

    Index id_readin = std::stoi(split[0]);
    if (id_readin != id) {
      throw std::runtime_error("Ids should be sorted and start at zero.");
    }
    id++;
    result.push_back(std::stod(split[1]));
  }
  return result;
}

std::vector<Rate_Engine::PairRates> VAverage::ReadRatefile(
    std::string filename) const {
  std::vector<Rate_Engine::PairRates> result;
  std::ifstream intt;
  intt.open(filename);
  if (!intt.is_open()) {
    throw std::runtime_error("File:" + filename + " could not be opened");
  }
  Index id = 0;
  std::string line;
  while (intt.good()) {

    tools::getline(intt, line);
    std::vector<std::string> split = tools::Tokenizer(line, " \t").ToVector();
    if (!split.size() || split[0] == "#" || split[0].substr(0, 1) == "#") {
      continue;
    }
    if (split.size() != 5) {
      throw std::runtime_error(
          "Row should only contain pairid,segid1,segid2 ,rate12,rate21");
    }

    Index id_readin = std::stoi(split[0]);
    if (id_readin != id) {
      throw std::runtime_error("Ids should be sorted and start at zero.");
    }
    id++;
    Rate_Engine::PairRates pair;
    pair.rate12 = std::stod(split[3]);
    pair.rate21 = std::stod(split[4]);
    result.push_back(pair);
  }
  return result;
}

bool VAverage::Evaluate(Topology& top) {
  std::cout << std::endl
            << "... ... Computing velocity average for all sites\n";
  std::cout << "Reading in site occupations from " << occfile_ << std::endl;
  std::vector<double> occ = ReadOccfile(occfile_);
  if (top.Segments().size() != occ.size()) {
    throw std::runtime_error(
        "Number of occupations is" + std::to_string(occ.size()) +
        " Topology has size:" + std::to_string(top.Segments().size()));
  }
  std::cout << "Reading in rates from " << ratefile_ << std::endl;
  std::vector<Rate_Engine::PairRates> rates = ReadRatefile(ratefile_);
  if (top.NBList().size() != Index(rates.size())) {
    throw std::runtime_error(
        "Number of pairs in file is" + std::to_string(rates.size()) +
        " Topology has size:" + std::to_string(top.NBList().size()));
  }
  std::vector<Eigen::Vector3d> velocities =
      std::vector<Eigen::Vector3d>(occ.size(), Eigen::Vector3d::Zero());
  for (const QMPair* pair : top.NBList()) {
    Index id1 = pair->Seg1()->getId();
    Index id2 = pair->Seg1()->getId();
    Index pairid = pair->getId();
    double p1 = occ[id1];
    double p2 = occ[id2];
    double w12 = rates[pairid].rate12;
    double w21 = rates[pairid].rate21;
    const Eigen::Vector3d& dR12 =
        pair->R();  // points from seg. (1) to seg. (2)
    // Split current symmetrically among sites (1) and (2)
    Eigen::Vector3d v_12_21 = 0.5 * (p1 * w12 - p2 * w21) * dR12;
    velocities[id1] += v_12_21;
    velocities[id2] += v_12_21;
  }

  std::ofstream ofs;
  ofs.open(outputfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + outputfile_);
  }
  ofs << "# 1   |   2     |  3 4 5       |   6  7  8         |" << std::endl;
  ofs << "# ----|---------|--------------|-------------------|" << std::endl;
  ofs << "# ID  |   NAME  |  X Y Z [nm]  |   VX VY VZ [nm/s]  |" << std::endl;
  for (const Segment& seg : top.Segments()) {
    const Eigen::Vector3d r = seg.getPos() * tools::conv::bohr2nm;
    const Eigen::Vector3d v = velocities[seg.getId()] * tools::conv::bohr2nm;
    ofs << (boost::format("%1$4d %2$-10s %3$+1.7e %4$+1.7e "
                          "%5$+1.7e %6$+1.7e %7$+1.7e %8$+1.7e") %
            seg.getId() % seg.getType() % r.x() % r.y() % r.z() % v.x() %
            v.y() % v.z())
        << std::endl;
  }
  std::cout << "Writing velocities to " << outputfile_ << std::endl;
  return true;
}

}  // namespace xtp
}  // namespace votca
