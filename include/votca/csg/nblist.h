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

#ifndef _VOTCA_CSG_NBLIST_H
#define _VOTCA_CSG_NBLIST_H

#include "beadlist.h"
#include "beadpair.h"
#include "exclusionlist.h"
#include "pairlist.h"

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;
/**
 * \brief Neighbour list class
 *
 * Implements a simple N^2 neighbour search and stores neigbourlist with pair
 * structure. User defined criteria can be added by SetMatchFunction. To only
 * get every pair listed once, the SetMatchFunction can be used and always
 * return that the pair is not stored.
 *
 */
class NBList : public PairList<Bead *, BeadPair> {
 public:
  NBList();
  virtual ~NBList();

  /// Generate the neighbour list based on two bead lists (e.g. bead types)
  virtual void Generate(BeadList &list1, BeadList &list2,
                        bool do_exclusions = true);
  /// Generate the neighbour list based on a single bead list
  virtual void Generate(BeadList &list, bool do_exclusions = true) {
    Generate(list, list, do_exclusions);
  }

  /// set the cutoff for the neighbour search
  void setCutoff(double cutoff) { _cutoff = cutoff; }
  /// get the cutoff for the neighbour search
  double getCutoff() { return _cutoff; }

  /**
   *  \brief match function for class member functions
   *
   * SetMatchFunction can be used to specify additional criteria, weather two
   * beads are added to the list of pairs or not. The function gets the two
   * two beads and the distance vector as argument. If a pair should be added,
   * the function should return true, otherwise false.
   *
   * This function can also be used, in a situation where each pair needs only
   * to be processed once, but the total number of pairs is to big to be stored
   * in memory, e.g. to calculate rdf for huge systems. In this case, set a
   * match function which always returns false (->no pair is added), and do
   * the processing in the match function.
   */
  template <typename T>
  void SetMatchFunction(T *object,
                        bool (T::*fkt)(Bead *, Bead *, const Eigen::Vector3d &,
                                       const double dist));

  /// \brief match function for static member functions or plain functions
  void SetMatchFunction(bool (*fkt)(Bead *, Bead *, const Eigen::Vector3d &,
                                    const double dist));

  /// standard match function
  static bool match_always(Bead *b1, Bead *b2, const Eigen::Vector3d &r,
                           const double dist) {
    return true;
  }

  /// function to use a user defined pair type
  template <typename pair_type>
  void setPairType();

 protected:
  /// cutoff
  double _cutoff;
  /// take into account exclusions from topolgoy
  bool _do_exclusions;

  /// policy function to create new bead types
  template <typename pair_type>
  static BeadPair *beadpair_create_policy(Bead *bead1, Bead *bead2,
                                          const Eigen::Vector3d &r) {
    return dynamic_cast<BeadPair *>(new pair_type(bead1, bead2, r));
  }

  typedef BeadPair *(*pair_creator_t)(Bead *bead1, Bead *bead2,
                                      const Eigen::Vector3d &r);
  /// the current bead pair creator function
  pair_creator_t _pair_creator;

 protected:
  /// Functor for match function to be able to set member and non-member
  /// functions
  class Functor {
   public:
    Functor() {}
    virtual bool operator()(Bead *, Bead *, const Eigen::Vector3d &,
                            const double dist) = 0;
    virtual ~Functor(){};
  };

  /// Functor for member functions
  template <typename T>
  class FunctorMember : public Functor {
   public:
    typedef bool (T::*fkt_t)(Bead *, Bead *, const Eigen::Vector3d &,
                             const double dist);

    FunctorMember(T *cls, fkt_t fkt) : _cls(cls), _fkt(fkt) {}

    bool operator()(Bead *b1, Bead *b2, const Eigen::Vector3d &r,
                    const double dist) {
      return (_cls->*_fkt)(b1, b2, r, dist);
    }

   private:
    T *   _cls;
    fkt_t _fkt;
  };

  /// Functor for non-member functions
  class FunctorNonMember : public Functor {
   public:
    typedef bool (*fkt_t)(Bead *, Bead *, const Eigen::Vector3d &,
                          const double dist);
    FunctorNonMember(fkt_t fkt) : _fkt(fkt) {}

    bool operator()(Bead *b1, Bead *b2, const Eigen::Vector3d &r,
                    const double dist) {
      return (*_fkt)(b1, b2, r, dist);
    }

   private:
    fkt_t _fkt;
  };

  Functor *_match_function;
};

template <typename pair_type>
void NBList::setPairType() {
  _pair_creator = NBList::beadpair_create_policy<pair_type>;
}

template <typename T>
inline void NBList::SetMatchFunction(T *object,
                                     bool (T::*fkt)(Bead *, Bead *,
                                                    const Eigen::Vector3d &,
                                                    const double)) {
  if (_match_function) delete _match_function;
  _match_function = dynamic_cast<Functor *>(new FunctorMember<T>(object, fkt));
}

inline void NBList::SetMatchFunction(bool (*fkt)(Bead *, Bead *,
                                                 const Eigen::Vector3d &,
                                                 const double)) {
  if (_match_function) delete _match_function;
  _match_function = dynamic_cast<Functor *>(new FunctorNonMember(fkt));
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_NBLIST_H */
