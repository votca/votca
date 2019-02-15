/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef __STATESERVER_H
#define __STATESERVER_H

#include <boost/format.hpp>
#include <votca/ctp/qmcalculator.h>

namespace votca {
namespace xtp {

class StateServer : public ctp::QMCalculator {
 public:
  StateServer(){};
  ~StateServer(){};

  string Identify() { return "stateserver"; }

  void Initialize(Property *options);
  bool EvaluateFrame(ctp::Topology *top);

  void DownloadTopology(FILE *out, ctp::Topology *top);
  void DownloadSegments(FILE *out, ctp::Topology *top);
  void DownloadPairs(FILE *out, ctp::Topology *top);
  void DownloadEList(FILE *out, ctp::Topology *top);
  void DownloadIList(FILE *out, ctp::Topology *top);
  void DownloadCoords(FILE *out, ctp::Topology *top){};

  void WriteXQM(FILE *out, ctp::Topology *top);
  void WriteEMP(FILE *out, ctp::Topology *top);
  void WriteEXCITED(FILE *out, ctp::Topology *top);
  void WriteULM(ctp::Topology *top);

 private:
  string _outfile;
  string _pdbfile;
  vector<string> _keys;

  string _xmp_alloc_file;
  string _emp_alloc_file;
};

void StateServer::Initialize(Property *opt) {

  // update options with the VOTCASHARE defaults
  UpdateWithDefaults(opt, "xtp");

  string tag = "options." + Identify();

  // Tabular output
  if (opt->exists(tag + ".out")) {
    _outfile = opt->get(tag + ".out").as<string>();
  } else {
    _outfile = "stateinfo.out";
  }
  // PDB output
  if (opt->exists(tag + ".pdb")) {
    _pdbfile = opt->get(tag + ".pdb").as<string>();
  } else {
    _pdbfile = "system.pdb";
  }

  string keys = "";
  if (opt->exists(tag + ".keys")) {
    keys = opt->get(tag + ".keys").as<string>();
  }

  Tokenizer tok_keys(keys, " ");
  tok_keys.ToVector(_keys);
}

bool StateServer::EvaluateFrame(ctp::Topology *top) {

  // ++++++++++++++++++++++++ //
  // Topology - Sites - Pairs //
  // ++++++++++++++++++++++++ //

  vector<string>::iterator key;
  bool writeTrajectory = false;

  string outfile =
      boost::lexical_cast<string>(top->getDatabaseId()) + "_" + _outfile;

  FILE *out = NULL;
  out = fopen(outfile.c_str(), "w");
  for (key = _keys.begin(); key < _keys.end(); key++) {

    cout << endl << "... ... Write " << flush;

    if (*key == "topology") {

      cout << "topology, ";

      fprintf(out, "           # +++++++++++++++++++ # \n");
      fprintf(out, "           # MD|QM Topology Info # \n");
      fprintf(out, "           # +++++++++++++++++++ # \n\n");

      DownloadTopology(out, top);
    } else if (*key == "sites") {

      cout << "sites, ";

      fprintf(out, "           # ++++++++++++++++++++ # \n");
      fprintf(out, "           # Segments in Database # \n");
      fprintf(out, "           # ++++++++++++++++++++ # \n\n");

      DownloadSegments(out, top);
    }

    else if (*key == "pairs") {

      cout << "pairs, ";

      fprintf(out, "           # +++++++++++++++++ # \n");
      fprintf(out, "           # Pairs in Database # \n");
      fprintf(out, "           # +++++++++++++++++ # \n\n");

      DownloadPairs(out, top);

    }

    else if (*key == "ilist") {

      cout << "integrals, ";

      fprintf(out, "           # +++++++++++++++++++++ # \n");
      fprintf(out, "           # Integrals in Database # \n");
      fprintf(out, "           # +++++++++++++++++++++ # \n\n");

      DownloadIList(out, top);

    }

    else if (*key == "elist") {

      cout << "energies, ";

      fprintf(out, "           # +++++++++++++++++++++++++ # \n");
      fprintf(out, "           # Site energies in Database # \n");
      fprintf(out, "           # +++++++++++++++++++++++++ # \n\n");

      DownloadEList(out, top);

    }

    else if (*key == "xqm") {

      cout << "XQMultipole Site Jobs, ";

      FILE *out_xqm;
      string job_file = "jobs.xml";
      out_xqm = fopen(job_file.c_str(), "w");

      WriteXQM(out_xqm, top);

      fclose(out_xqm);
    }

    else if (*key == "mps") {

      cout << "MPS background, ";

      FILE *out_emp;
      string emp_file = "mps.tab";
      out_emp = fopen(emp_file.c_str(), "w");

      WriteEMP(out_emp, top);

      fclose(out_emp);
    }

    else {
      cout << "ERROR (Invalid key " << *key << ") ";
    }

    cout << "done. " << flush;

    fprintf(out, "\n\n");
  }
  fclose(out);

  if (writeTrajectory) {

    // ++++++++++++++++++ //
    // MD, QM Coordinates //
    // ++++++++++++++++++ //
    // Segments
    string pdbfile = boost::lexical_cast<string>(top->getDatabaseId()) +
                     "_conjg_" + _pdbfile;

    out = fopen(pdbfile.c_str(), "w");
    top->WritePDB(out);
    fclose(out);

    // Fragments
    pdbfile = boost::lexical_cast<string>(top->getDatabaseId()) + "_rigid_" +
              _pdbfile;

    out = fopen(pdbfile.c_str(), "w");
    vector<ctp::Segment *>::iterator segIt;
    for (segIt = top->Segments().begin(); segIt < top->Segments().end();
         segIt++) {
      ctp::Segment *seg = *segIt;
      seg->WritePDB(out);
    }
    fclose(out);
  }

  return true;
}

void StateServer::DownloadTopology(FILE *out, ctp::Topology *top) {

  fprintf(out, "Topology Database ID %3d \n", top->getDatabaseId());
  fprintf(out,
          "  Periodic Box: %2.4f %2.4f %2.4f | %2.4f %2.4f %2.4f "
          "| %2.4f %2.4f %2.4f \n",
          top->getBox().get(0, 0), top->getBox().get(0, 1),
          top->getBox().get(0, 2), top->getBox().get(1, 0),
          top->getBox().get(1, 1), top->getBox().get(1, 2),
          top->getBox().get(2, 0), top->getBox().get(2, 1),
          top->getBox().get(2, 2));

  fprintf(out, "  Step number %7d \n", top->getStep());
  fprintf(out, "  Time          %2.3f \n", top->getTime());

  int N_mol = top->Molecules().size();
  int N_seg = top->Segments().size();
  int N_atm = top->Atoms().size();
  int N_nbs = top->NBList().size();

  fprintf(out, "  # Molecules %7d \n", N_mol);
  fprintf(out, "  # Segments  %7d \n", N_seg);
  fprintf(out, "  # Atoms     %7d \n", N_atm);
  fprintf(out, "  # Pairs     %7d \n", N_nbs);
  return;
}

void StateServer::DownloadSegments(FILE *out, ctp::Topology *top) {

  vector<ctp::Segment *>::iterator segit;

  for (segit = top->Segments().begin(); segit < top->Segments().end();
       segit++) {

    ctp::Segment *seg = *segit;
    fprintf(out,
            "SiteID %5d %5s | xyz %8.3f %8.3f %8.3f "
            "| SiteE(intra) C %2.4f A %2.4f S %2.4f T %2.4f "
            "| Lambdas: NC %2.4f CN %2.4f NA %2.4f AN %2.4f "
            "| Lambdas_ex: NS %2.4f SN %2.4f NT %2.4f TN %2.4f "
            "| SiteE(pol+estat) C %2.4f N %2.4f A %2.4f S %2.4f T %2.4f    \n",
            seg->getId(), seg->getName().c_str(), seg->getPos().getX(),
            seg->getPos().getY(), seg->getPos().getZ(), seg->getU_cC_nN(1),
            seg->getU_cC_nN(-1), seg->getU_xX_nN(+2), seg->getU_xX_nN(+3),
            seg->getU_nC_nN(1), seg->getU_cN_cC(1), seg->getU_nC_nN(-1),
            seg->getU_cN_cC(-1), seg->getU_nX_nN(+2), seg->getU_xN_xX(+2),
            seg->getU_nX_nN(+3), seg->getU_xN_xX(+3), seg->getEMpoles(1),
            seg->getEMpoles(0), seg->getEMpoles(-1), seg->getEMpoles(2),
            seg->getEMpoles(3));
  }
  return;
}

void StateServer::DownloadPairs(FILE *out, ctp::Topology *top) {
  ctp::QMNBList::iterator nit;

  for (nit = top->NBList().begin(); nit != top->NBList().end(); nit++) {

    ctp::QMPair *pair = *nit;

    int ghost = (pair->HasGhost()) ? 1 : 0;

    fprintf(out,
            "ID %5d SEG1 %4d %1s SEG2 %4d %1s dR %2.4f PBC %1d "
            "dE(-1) %4.7f dE(+1) %4.7f "
            "L(-1) %1.4f L(+1) %1.4f J2(-1) %1.8e J2(+1) %1.8e "
            "R12(-1) %2.4e R12(+1) %2.4e R21(-1) %2.4e R21(-1) %2.4e "
            "dE(+2) %4.7f dE(+3) %4.7f "
            "L(+2) %1.4f L(+3) %1.4f J2(+2) %1.8e J2(+3) %1.8e "
            "R12(+2) %2.4e R12(+3) %2.4e R21(+2) %2.4e R21(+3) %2.4e\n",
            pair->getId(), pair->first->getId(), pair->first->getName().c_str(),
            pair->second->getId(), pair->second->getName().c_str(),
            pair->Dist(), ghost, pair->getdE12(-1), pair->getdE12(+1),
            pair->getLambdaO(-1), pair->getLambdaO(+1), pair->getJeff2(-1),
            pair->getJeff2(+1), pair->getRate12(-1), pair->getRate12(+1),
            pair->getRate21(-1), pair->getRate21(+1), pair->getdE12(+2),
            pair->getdE12(+3), pair->getLambdaO(+2), pair->getLambdaO(+3),
            pair->getJeff2(+2), pair->getJeff2(+3), pair->getRate12(+2),
            pair->getRate12(+3), pair->getRate21(+2), pair->getRate21(+3));
  }
  return;
}

void StateServer::DownloadIList(FILE *out, ctp::Topology *top) {
  ctp::QMNBList::iterator nit;
  for (nit = top->NBList().begin(); nit != top->NBList().end(); ++nit) {
    ctp::QMPair *pair = *nit;

    fprintf(out,
            "%5d %5d %5d e %4.7e h %4.7e s %4.7e t %4.7e edr %4.7f pbc %1i\n",
            pair->getId(), pair->Seg1()->getId(), pair->Seg2()->getId(),
            pair->getJeff2(-1), pair->getJeff2(+1), pair->getJeff2(+2),
            pair->getJeff2(+3), pair->Dist(), (pair->HasGhost()) ? 1 : 0);
  }
  return;
}

void StateServer::WriteEMP(FILE *out, ctp::Topology *top) {

  fprintf(out, "# ID   TYPE    _n.mps    _e.mps    _h.mps \n");

  vector<Segment *>::iterator sit;
  for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

    fprintf(out, "%4d %15s %-30s %-30s %-30s \n", (*sit)->getId(),
            (*sit)->getName().c_str(),
            ("MP_FILES/" + (*sit)->getName() + "_n.mps").c_str(),
            ("MP_FILES/" + (*sit)->getName() + "_e.mps").c_str(),
            ("MP_FILES/" + (*sit)->getName() + "_h.mps").c_str());
  }
  return;
}

void StateServer::WriteXQM(FILE *out, ctp::Topology *top) {

  fprintf(out, "<jobs>\n");

  //    QMNBList::iterator nit;
  //    for (nit = top->NBList().begin();
  //         nit != top->NBList().end();
  //         nit++) {
  //
  //        QMPair *qmpair = *nit;
  //
  //       string prefix = "pair_"+boost::lexical_cast<string>(qmpair->getId())
  //                   +"_"+boost::lexical_cast<string>(qmpair->Seg1()->getId())
  //                   +"_"+boost::lexical_cast<string>(qmpair->Seg2()->getId());
  //       string tag = "tag_xxx";
  //
  //        fprintf(out, "%5d %5s   %5d    %4d %5s %-30s   %4d %5s %-30s \n",
  //                     qmpair->getId(),
  //                     "tag_xxx",
  //                     qmpair->getId(),
  //                     qmpair->first->getId(),
  //                     qmpair->first->getName().c_str(),
  //                     (prefix+"_1.mps").c_str(),
  //                     qmpair->second->getId(),
  //                     qmpair->second->getName().c_str(),
  //                     (prefix+"_2.mps").c_str());
  //    }

  int jobId = 0;
  vector<ctp::Segment *>::iterator sit;
  for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

    ctp::Segment *seg = *sit;
    int segId = seg->getId();
    string segName = seg->getName();

    string stateStr = "n";
    string jobTag = boost::lexical_cast<string>(segId) + "_" + stateStr;
    string mpsFile = "MP_FILES/" + segName + "_" + stateStr + ".mps";
    ++jobId;
    string jobStr = (boost::format("<job>\n"
                                   "\t<id>%1$d</id>\n"
                                   "\t<tag>%2$s</tag>\n"
                                   "\t<input>%3$d:%4$s:%5$s</input>\n"
                                   "\t<status>AVAILABLE</status>\n"
                                   "</job>\n") %
                     jobId % jobTag % segId % segName % mpsFile)
                        .str();
    fprintf(out, "%s", jobStr.c_str());

    stateStr = "e";
    jobTag = boost::lexical_cast<string>(segId) + "_" + stateStr;
    mpsFile = "MP_FILES/" + segName + "_" + stateStr + ".mps";
    ++jobId;
    jobStr = (boost::format("<job>\n"
                            "\t<id>%1$d</id>\n"
                            "\t<tag>%2$s</tag>\n"
                            "\t<input>%3$d:%4$s:%5$s</input>\n"
                            "\t<status>AVAILABLE</status>\n"
                            "</job>\n") %
              jobId % jobTag % segId % segName % mpsFile)
                 .str();
    fprintf(out, "%s", jobStr.c_str());

    stateStr = "h";
    jobTag = boost::lexical_cast<string>(segId) + "_" + stateStr;
    mpsFile = "MP_FILES/" + segName + "_" + stateStr + ".mps";
    ++jobId;
    jobStr = (boost::format("<job>\n"
                            "\t<id>%1$d</id>\n"
                            "\t<tag>%2$s</tag>\n"
                            "\t<input>%3$d:%4$s:%5$s</input>\n"
                            "\t<status>AVAILABLE</status>\n"
                            "</job>\n") %
              jobId % jobTag % segId % segName % mpsFile)
                 .str();
    fprintf(out, "%s", jobStr.c_str());
  }

  fprintf(out, "</jobs>\n");
  return;
}

}  // namespace xtp
}  // namespace votca

#endif
