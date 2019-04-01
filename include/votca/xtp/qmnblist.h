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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef _VOTCA_XTP_QMNBList_H
#define _VOTCA_XTP_QMNBList_H

#include <list>
#include <string>
#include <vector>

#include <votca/csg/pairlist.h>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/qmpair.h>

namespace CSG = votca::csg;

namespace votca {
namespace xtp {

class Topology;

class QMNBList : public CSG::PairList<Segment *, QMPair> {
 public:
  /**
   * \brief Container for records of type Donor-Bridge1-Bridge2-...-Acceptor
   *
   * Every SuperExchangeType record contains a pair donor, acceptor
   * and a list of bridges (all specified by segment types, i.e. strings)
   *
   */
  class SuperExchangeType {
   public:
    // Initializes the object from a [Donor Bridge1 Bridge2 ... Acceptor] string
    SuperExchangeType(std::string initString) {

      votca::tools::Tokenizer tok(initString, " ");
      std::vector<std::string> names;
      tok.ToVector(names);

      if (names.size() < 3) {
        std::cout << "ERROR: Faulty superexchange definition: "
                  << "Need at least three segment names (DONOR BRIDGES "
                     "ACCEPTOR separated by a space"
                  << std::endl;
        throw std::runtime_error("Error in options file.");
      }

      // fill up the donor-bride-acceptor structure

      donor = names.front();
      acceptor = names.back();

      for (std::vector<std::string>::iterator it = ++names.begin();
           it != --names.end(); it++) {
        bridges.push_back(*it);
      }
    }

    bool isOfBridge(std::string segment_type) {
      std::list<std::string>::iterator findIter =
          std::find(bridges.begin(), bridges.end(), segment_type);
      return findIter != bridges.end();
    };

    bool isOfDonorAcceptor(std::string segment_type) {
      return segment_type == donor || segment_type == acceptor;
    }

    std::string asString() {
      std::string ts;
      ts += donor;
      for (std::list<std::string>::iterator si = bridges.begin();
           si != bridges.end(); si++)
        ts = ts + " " + *si;
      ts += " " + acceptor;
      return ts;
    }

   private:
    std::string donor;
    std::string acceptor;
    std::list<std::string> bridges;
  };

  QMNBList() : _top(NULL), _cutoff(0){};
  QMNBList(Topology *top) : _top(top), _cutoff(0){};
  ~QMNBList() {
    CSG::PairList<Segment *, QMPair>::Cleanup();
    // cleanup the list of superexchange pairs
    for (std::list<SuperExchangeType *>::iterator it = _superexchange.begin();
         it != _superexchange.end(); it++) {
      delete *it;
    }
  }

  /**
   * \brief Adds SuperExchange pairs to the neighbor list
   *
   * SuperExchange pairs are those pairs which have one or more bridging
   * molecules specified in the input file via a record [Donor Bridge1 Bridge2
   * ... Acceptor] Every pair gets a flag, identifying if a QMPair is of type
   * Hopping, SuperExchange, SuperExchangeAndHopping. This is stored to the
   * state file The BRIDGED pairs are stored but BRIDGING pairs have to be
   * regenerated every time we need them (edft job writer, idft job writer and
   * importer)
   *
   */
  void GenerateSuperExchange();

  /**
   * @param type Adds a SuperExchangeType based on this string (Donor Bridge1
   * Bridge2 ... Acceptor)
   */
  void AddSuperExchangeType(std::string type) {
    _superexchange.push_back(new SuperExchangeType(type));
  }

  void setSuperExchangeTypes(std::list<SuperExchangeType *> types) {
    _superexchange = types;
  }

  const std::list<SuperExchangeType *> &getSuperExchangeTypes() const {
    return _superexchange;
  }

  void setCutoff(double cutoff) { _cutoff = cutoff; }
  double getCutoff() { return _cutoff; }

  QMPair *Add(Segment *seg1, Segment *seg2, bool safe = true);

  void PrintInfo(std::FILE *out);

  void AddQMNBlist(QMNBList &temp);

 protected:
  Topology *_top;
  double _cutoff;
  std::list<SuperExchangeType *> _superexchange;
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_QMNBLIST_H */
