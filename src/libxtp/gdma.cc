/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/ctp/logger.h>
#include <votca/xtp/gdma.h>

namespace votca {
namespace xtp {

using namespace std;

// initialize the GDMA object, set parameters
// for use of external code -> Rank

void GDMA::Initialize(tools::Property& options) {

  string key = "gdma";
  _chkFile = options.ifExistsReturnElseReturnDefault<string>(key + ".chk",
                                                             "system.chk");
  _executable = options.ifExistsReturnElseReturnDefault<string>(
      key + ".executable", "gdma");
  _density =
      options.ifExistsReturnElseReturnDefault<string>(key + ".density", "SCF");
  _limit = options.ifExistsReturnElseReturnDefault<int>(
      key + ".multipoles.limit", 2);

  if (_limit > 2) {
    throw std::runtime_error("Tried to use GDMA with Rank > 2. Not supported!");
  }

  _switch = options.ifExistsReturnElseReturnDefault<double>(
      key + ".multipoles.switch", 4.0);  // corresponds to GDMA defaul
  _outFile = options.ifExistsReturnElseReturnDefault<string>(key + ".output",
                                                             "gdma.out");
}

// write an input file for the external GDMA code by A. Stone
void GDMA::WriteInputFile() {
  ofstream gdma_inputfile;
  string gdma_inputfile_name_full = _runFolder + "/gdma.in";
  gdma_inputfile.open(gdma_inputfile_name_full.c_str());
  gdma_inputfile << "Title \"Multipole Fit for QMMM\"" << endl;
  gdma_inputfile << "File system.fchk Density " << _density << endl;
  gdma_inputfile << "Angstrom" << endl;
  gdma_inputfile << "Multipoles" << endl;
  gdma_inputfile << "Limit " << _limit << endl;
  gdma_inputfile << "Switch " << _switch << endl;
  gdma_inputfile << "Start" << endl;
  gdma_inputfile << "Finish" << endl;
  gdma_inputfile.close();
}

void GDMA::RunExternal() {

  // check if the input file exists
  string fullInput = _runFolder + "/gdma.in";
  if (!boost::filesystem::exists(fullInput)) {
    CTP_LOG(ctp::logINFO, *_log)
        << "GDMA input file has not been found!" << flush;
    throw runtime_error(" GDMA cannot be run! ");
  }

  // check if fchk exists
  string fullFChk = _runFolder + "/system.fchk";
  if (!boost::filesystem::exists(fullFChk)) {
    // try converting Chk to FChk
    string fullChk = _runFolder + "/" + _chkFile;
    if (boost::filesystem::exists(fullChk)) {
      // use formchk
      string command;
      command = "cd " + _runFolder + "; formchk " + _chkFile +
                " system.fchk > /dev/null";
      if (std::system(command.c_str())) {
        throw runtime_error("Command " + command + "failed");
      }
      // check again for fchk
      if (!boost::filesystem::exists(fullFChk)) {
        CTP_LOG(ctp::logINFO, *_log) << "Formatted Checkpoint file has not "
                                        "been found and cannot be created!"
                                     << flush;
        throw runtime_error(" GDMA cannot be run! ");
      }
    } else {
      CTP_LOG(ctp::logINFO, *_log) << "Formatted Checkpoint file has not been "
                                      "found and cannot be created!"
                                   << flush;
      throw runtime_error(" GDMA cannot be run! ");
    }
  }

  // now we seem ready to go
  string command;
  command =
      "cd " + _runFolder + "; " + _executable + " < gdma.in > " + _outFile;
  if (std::system(command.c_str())) {
    throw runtime_error("Command " + command + "failed");
  }
}

void GDMA::ParseOutputFile() {

  string gdma_output_name_full = _runFolder + "/gdma.out";
  std::ifstream gdma_output(gdma_output_name_full.c_str());
  std::string line;
  while (gdma_output) {

    getline(gdma_output, line);
    // if a line has an equality sign, must be energy
    std::string::size_type atom_pos = line.find("x =");
    if (atom_pos != std::string::npos) {
      // found an atom, read one more line
      getline(gdma_output, line);
      // determine rank
      std::vector<string> results;
      boost::trim(line);
      std::vector<double> Qs;  // temp vector for reading in

      boost::algorithm::split(results, line, boost::is_any_of("\t "),
                              boost::algorithm::token_compress_on);

      int rank = boost::lexical_cast<int>(results[3]);
      if (rank < 0 || rank > 2) {
        throw runtime_error(
            (boost::format(" Invalid GDMA rank %s!") % rank).str());
      }

      // atomic charge
      if (rank >= 0) {
        // getting charge
        getline(gdma_output, line);
        std::vector<string> results;
        boost::trim(line);
        boost::algorithm::split(results, line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
        double Q00 = boost::lexical_cast<double>(results.back());
        Qs.push_back(Q00);
      }

      // atomic dipole components
      if (rank >= 1) {

        // getting dipoles
        getline(gdma_output, line);
        std::vector<string> results;
        boost::trim(line);
        boost::algorithm::split(results, line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
        double Q10 = boost::lexical_cast<double>(results[5]);
        double Q11c = boost::lexical_cast<double>(results[8]);
        double Q11s = boost::lexical_cast<double>(results.back());
        Qs.push_back(Q10);
        Qs.push_back(Q11c);
        Qs.push_back(Q11s);
      }

      // atomic quadrupole components
      if (rank == 2) {

        getline(gdma_output, line);
        std::vector<string> results;
        boost::trim(line);
        boost::algorithm::split(results, line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
        double Q20 = boost::lexical_cast<double>(results[5]);
        double Q21c = boost::lexical_cast<double>(results[8]);
        double Q21s = boost::lexical_cast<double>(results.back());
        getline(gdma_output, line);
        boost::trim(line);
        boost::algorithm::split(results, line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
        double Q22c = boost::lexical_cast<double>(results[2]);
        double Q22s = boost::lexical_cast<double>(results[5]);

        Qs.push_back(Q20);
        Qs.push_back(Q21c);
        Qs.push_back(Q21s);
        Qs.push_back(Q22c);
        Qs.push_back(Q22s);
      }
      _multipoles.push_back(Qs);

    }  // atom
    std::string::size_type break_pos =
        line.find("Total multipoles referred to origin at");
    if (break_pos != std::string::npos) break;
  }  // gdma_output
}

}  // namespace xtp
}  // namespace votca
