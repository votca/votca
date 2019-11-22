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

#include "neighborlist.h"

using namespace std;

namespace votca {
namespace xtp {

void Neighborlist::Initialize(tools::Property& options) {

  // update options with the VOTCASHARE defaults
  UpdateWithDefaults(options, "xtp");
  std::string key = "options." + Identify();

  std::vector<tools::Property*> segs = options.Select(key + ".segments");

  for (tools::Property* segprop : segs) {
    std::string types = segprop->get("type").as<std::string>();
    double cutoff = segprop->get("cutoff").as<double>() * tools::conv::nm2bohr;

    tools::Tokenizer tok(types, " ");
    std::vector<std::string> names;
    tok.ToVector(names);

    if (names.size() != 2) {
      throw std::runtime_error(
          "ERROR: Faulty pair definition for cut-off's: Need two segment names "
          "separated by a space");
    }
    _cutoffs[names[0]][names[1]] = cutoff * tools::conv::nm2bohr;
    _cutoffs[names[1]][names[0]] = cutoff * tools::conv::nm2bohr;
    if (std::find(_included_segments.begin(), _included_segments.end(),
                  names[0]) == _included_segments.end()) {
      _included_segments.push_back(names[0]);
    }
    if (std::find(_included_segments.begin(), _included_segments.end(),
                  names[1]) == _included_segments.end()) {
      _included_segments.push_back(names[1]);
    }
  }

  if (options.exists(key + ".constant")) {
    _useConstantCutoff = true;
    _constantCutoff =
        options.get(key + ".constant").as<double>() * tools::conv::nm2bohr;
  } else {
    _useConstantCutoff = false;
  }
  if (options.exists(key + ".exciton_cutoff")) {
    _useExcitonCutoff = true;
    _excitonqmCutoff = options.get(key + ".exciton_cutoff").as<double>() *
                       tools::conv::nm2bohr;
  } else {
    _useExcitonCutoff = false;
  }
}

Index Neighborlist::DetClassicalPairs(Topology& top) {
  Index classical_pairs = 0;
#pragma omp parallel for
  for (Index i = 0; i < top.NBList().size(); i++) {
    const Segment* seg1 = top.NBList()[i]->Seg1();
    const Segment* seg2 = top.NBList()[i]->Seg2();
    if (top.GetShortestDist(*seg1, *seg2) > _excitonqmCutoff) {
      top.NBList()[i]->setType(QMPair::Excitoncl);
#pragma omp critical
      { classical_pairs++; }
    } else {
      top.NBList()[i]->setType(QMPair::Hopping);
    }
  }  // Type 3 Exciton_classical approx
  return classical_pairs;
}

bool Neighborlist::EvaluateFrame(Topology& top) {
  OPENMP::setMaxThreads(_nThreads);
  std::cout << " Using " << OPENMP::getMaxThreads() << " threads" << std::flush;

  if (tools::globals::verbose) {
    std::cout << std::endl << "... ..." << std::flush;
  }

  double min = top.getBox().diagonal().minCoeff();

  std::vector<Segment*> segs;
  for (Segment& seg : top.Segments()) {
    if (_useConstantCutoff ||
        std::find(_included_segments.begin(), _included_segments.end(),
                  seg.getType()) != _included_segments.end()) {
      segs.push_back(&seg);
      seg.getApproxSize();
    }
  }
  std::cout << std::endl;
  std::cout << "Evaluating " << segs.size() << " segments for neighborlist. ";
  if ((top.Segments().size() - segs.size()) != 0) {
    std::cout << top.Segments().size() - segs.size()
              << " segments are not taken into account as specified"
              << std::endl;
  } else {
    std::cout << std::endl;
  }
  if (!_useConstantCutoff) {
    std::cout << "The following segments are used in the neigborlist creation"
              << std::endl;
    std::cout << "\t" << std::flush;
    for (const std::string& st : _included_segments) {
      std::cout << " " << st << std::flush;
    }
    std::cout << std::endl;
  }

  std::cout << "\r ... ... Evaluating " << std::flush;
  std::vector<std::string> skippedpairs;

  top.NBList().Cleanup();

  boost::progress_display progress(segs.size());
  // cache approx sizes
  std::vector<double> approxsize = std::vector<double>(segs.size(), 0.0);
#pragma omp parallel for
  for (Index i = 0; i < Index(segs.size()); i++) {
    approxsize[i] = segs[i]->getApproxSize();
  }
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < Index(segs.size()); i++) {
    Segment* seg1 = segs[i];
    double cutoff = _constantCutoff;
    for (Index j = i + 1; j < Index(segs.size()); j++) {
      Segment* seg2 = segs[j];
      if (!_useConstantCutoff) {
        try {
          cutoff = _cutoffs.at(seg1->getType()).at(seg2->getType());
        } catch (const std::exception&) {
          std::string pairstring = seg1->getType() + "/" + seg2->getType();
          if (std::find(skippedpairs.begin(), skippedpairs.end(), pairstring) ==
              skippedpairs.end()) {
#pragma omp critical
            { skippedpairs.push_back(pairstring); }
          }
          continue;
        }
      }

      if (cutoff > 0.5 * min) {
        throw std::runtime_error(
            (boost::format("Cutoff is larger than half the box size. Maximum "
                           "allowed cutoff is %1$1.1f") %
             (tools::conv::nm2bohr * 0.5 * min))
                .str());
      }
      double cutoff2 = cutoff * cutoff;
      Eigen::Vector3d segdistance =
          top.PbShortestConnect(seg1->getPos(), seg2->getPos());
      double segdistance2 = segdistance.squaredNorm();
      double outside = cutoff + approxsize[i] + approxsize[j];

      if (segdistance2 < cutoff2) {
#pragma omp critical
        { top.NBList().Add(*seg1, *seg2, segdistance); }

      } else if (segdistance2 > (outside * outside)) {
        continue;
      } else {
        double R = top.GetShortestDist(*seg1, *seg2);
        if ((R * R) < cutoff2) {
#pragma omp critical
          { top.NBList().Add(*seg1, *seg2, segdistance); }
        }
      }
    } /* exit loop seg2 */
#pragma omp critical
    { ++progress; }
  } /* exit loop seg1 */

  if (skippedpairs.size() > 0) {
    std::cout << "WARNING: No cut-off specified for segment pairs of type "
              << std::endl;
    for (const std::string& st : skippedpairs) {
      std::cout << st << std::endl;
    }
    std::cout << "pairs were skipped" << std::endl;
  }

  std::cout << std::endl
            << " ... ... Created " << top.NBList().size() << " direct pairs.";
  if (_useExcitonCutoff) {
    std::cout << std::endl
              << " ... ... Determining classical pairs " << std::endl;
    Index classical_pairs = DetClassicalPairs(top);
    std::cout << " ... ... Found " << classical_pairs << " classical pairs "
              << std::endl;
  }

  // sort qmpairs by seg1id and then by seg2id then reindex the pair id
  // according to that.
  top.NBList().sortAndReindex([](QMPair* a, QMPair* b) {
    if (a->Seg1()->getId() != b->Seg1()->getId()) {
      return a->Seg1()->getId() < b->Seg1()->getId();
    }
    return a->Seg2()->getId() < b->Seg2()->getId();
  });

  return true;
}
}  // namespace xtp
}  // namespace votca
