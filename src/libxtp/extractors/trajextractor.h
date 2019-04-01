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

#ifndef VOTCA_XTP_TRAJEXTRACTOR_H
#define VOTCA_XTP_TRAJEXTRACTOR_H

#include <cstdlib>
#include <votca/ctp/qmcalculator.h>

namespace votca {
namespace xtp {

class TrajExtractor : public ctp::QMCalculator {

 public:
  TrajExtractor(){};
  ~TrajExtractor(){};

  std::string Identify() { return "extract.trajectory"; }

  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
  std::string _outPDBmd;
  std::string _outPDBqm;
};

void TrajExtractor::Initialize(tools::Property *options) {

  _outPDBmd = Identify() + "_md.pdb";
  _outPDBqm = Identify() + "_qm.pdb";
}

bool TrajExtractor::EvaluateFrame(ctp::Topology *top) {

  // Rigidify std::system (if possible)
  if (!top->Rigidify()) return 0;

  // Print coordinates
  FILE *outPDBmd = NULL;
  FILE *outPDBqm = NULL;

  outPDBmd = fopen(_outPDBmd.c_str(), "w");
  outPDBqm = fopen(_outPDBqm.c_str(), "w");

  std::fprintf(outPDBmd, "TITLE     VOT CAtastrophic title \n");
  std::fprintf(outPDBmd, "Model %8d \n", 1);

  std::fprintf(outPDBqm, "TITLE     VOT CAtastrophic title \n");
  std::fprintf(outPDBqm, "Model %8d \n", 1);

  std::vector<ctp::Segment *>::iterator sit;
  for (sit = top->Segments().begin(); sit < top->Segments().end(); sit++) {

    (*sit)->WritePDB(outPDBmd, "Atoms", "MD");
    (*sit)->WritePDB(outPDBqm, "Atoms", "QM");
  }

  std::fprintf(outPDBmd, "TER\nENDMDL\n");
  std::fprintf(outPDBqm, "TER\nENDMDL\n");

  std::fclose(outPDBmd);
  std::fclose(outPDBqm);

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TRAJEXTRACTOR_H
