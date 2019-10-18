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

#ifndef _VOTCA_CSG_NBLIST_3BODY_H
#define _VOTCA_CSG_NBLIST_3BODY_H

#include "beadlist.h"
#include "beadtriple.h"
#include "exclusionlist.h"
#include "triplelist.h"

namespace votca {
namespace csg {

/**
 * \brief Neighbour list class for 3 body interactions
 *
 * Implements a simple N^3 neighbour search and stores 3body neigbourlist with
 * triple structure.
 *
 *
 * User defined criteria can be added by SetMatchFunction. To only
 * get every pair listed once, the SetMatchFunction can be used and always
 * return that the pair is not stored.
 *
 */
class NBList_3Body : public TripleList<Bead *, BeadTriple> {
 public:
  NBList_3Body();
  ~NBList_3Body() override;

  /// Generate the 3body neighbour list based on three bead lists (e.g. bead
  /// types)
  virtual void Generate(BeadList &list1, BeadList &list2, BeadList &list3,
                        bool do_exclusions = true);
  /// Generate the 3body neighbour list based on two bead lists (e.g. bead
  /// types) Experimental: here the second and the third bead list are the same
  /// (nblist will have the structure bead type 1, bead type 2, bead type 2)
  virtual void Generate(BeadList &list1, BeadList &list2,
                        bool do_exclusions = true) {
    Generate(list1, list2, list2, do_exclusions);
  };
  /// Generate the 3body neighbour list based on a single bead list (nblist will
  /// have the structure bead type 1, bead type 1, bead type 1)
  virtual void Generate(BeadList &list, bool do_exclusions = true) {
    Generate(list, list, list, do_exclusions);
  }

  /// set the cutoff for the neighbour search
  /// to do: at the moment use only one single cutoff value
  void setCutoff(const double cutoff) { _cutoff = cutoff; }
  /// get the cutoff for the neighbour search
  double getCutoff() { return _cutoff; }

  /**
   *  \brief match function for class member functions
   *
   * SetMatchFunction can be used to specify additional criteria, weather three
   * beads are added to the list of triples or not. The function gets the
   * three beads and the distance vectors between the beads as argument.
   * If a triple should be added, the function should return true, otherwise
   * false.
   *
   * This function can also be used, in a situation where each triple needs only
   * to be processed once, but the total number of triples is to big to be
   * stored in memory, e.g. to calculate rdf for huge systems. In this case, set
   * a match function which always returns false (->no triple is added), and do
   * the processing in the match function.
   */
  template <typename T>
  void SetMatchFunction(
      T *object, bool (T::*fkt)(Bead *, Bead *, Bead *, const Eigen::Vector3d &,
                                const Eigen::Vector3d &,
                                const Eigen::Vector3d &, const double dist12,
                                const double dist13, const double dist23));

  /// \brief match function for static member functions or plain functions
  void SetMatchFunction(bool (*fkt)(
      Bead *, Bead *, Bead *, const Eigen::Vector3d &, const Eigen::Vector3d &,
      const Eigen::Vector3d &, const double dist12, const double dist13,
      const double dist23));

  /// standard match function
  static bool match_always(Bead *b1, Bead *b2, Bead *b3,
                           const Eigen::Vector3d &r12,
                           const Eigen::Vector3d &r13,
                           const Eigen::Vector3d &r23, const double dist12,
                           const double dist13, const double dist23) {
    return true;
  }

  /// function to use a user defined triple type
  template <typename triple_type>
  void setTripleType();

 protected:
  /// cutoff (at the moment use only one cutoff value)
  double _cutoff;
  /// take into account exclusions from topolgoy
  bool _do_exclusions;

  /// policy function to create new bead types
  template <typename triple_type>
  static BeadTriple *beadtriple_create_policy(Bead *bead1, Bead *bead2,
                                              Bead *bead3,
                                              const Eigen::Vector3d &r12,
                                              const Eigen::Vector3d &r13,
                                              const Eigen::Vector3d &r23) {
    return dynamic_cast<BeadTriple *>(
        new triple_type(bead1, bead2, bead3, r12, r13, r23));
  }

  using triple_creator_t = BeadTriple *(*)(Bead *, Bead *, Bead *,
                                           const Eigen::Vector3d &,
                                           const Eigen::Vector3d &,
                                           const Eigen::Vector3d &);
  /// the current bead pair creator function
  triple_creator_t _triple_creator;

 protected:
  /// Functor for match function to be able to set member and non-member
  /// functions
  class Functor {
   public:
    Functor() = default;
    virtual bool operator()(Bead *, Bead *, Bead *, const Eigen::Vector3d &,
                            const Eigen::Vector3d &, const Eigen::Vector3d &,
                            const double dist12, const double dist13,
                            const double dist23) = 0;
    virtual ~Functor() = default;
    ;
  };

  /// Functor for member functions
  template <typename T>
  class FunctorMember : public Functor {
   public:
    using fkt_t = bool (T::*)(Bead *, Bead *, Bead *, const Eigen::Vector3d &,
                              const Eigen::Vector3d &, const Eigen::Vector3d &,
                              const double, const double, const double);

    FunctorMember(T *cls, fkt_t fkt) : _cls(cls), _fkt(fkt) {}

    bool operator()(Bead *b1, Bead *b2, Bead *b3, const Eigen::Vector3d &r12,
                    const Eigen::Vector3d &r13, const Eigen::Vector3d &r23,
                    const double dist12, const double dist13,
                    const double dist23) override {
      return (_cls->*_fkt)(b1, b2, b3, r12, r13, r23, dist12, dist13, dist23);
    }

   private:
    T *_cls;
    fkt_t _fkt;
  };

  /// Functor for non-member functions
  class FunctorNonMember : public Functor {
   public:
    using fkt_t = bool (*)(Bead *, Bead *, Bead *, const Eigen::Vector3d &,
                           const Eigen::Vector3d &, const Eigen::Vector3d &,
                           const double, const double, const double);
    FunctorNonMember(fkt_t fkt) : _fkt(fkt) {}

    bool operator()(Bead *b1, Bead *b2, Bead *b3, const Eigen::Vector3d &r12,
                    const Eigen::Vector3d &r13, const Eigen::Vector3d &r23,
                    const double dist12, const double dist13,
                    const double dist23) override {
      return (*_fkt)(b1, b2, b3, r12, r13, r23, dist12, dist13, dist23);
    }

   private:
    fkt_t _fkt;
  };

  Functor *_match_function;
};

template <typename triple_type>
void NBList_3Body::setTripleType() {
  _triple_creator = NBList_3Body::beadtriple_create_policy<triple_type>;
}

template <typename T>
inline void NBList_3Body::SetMatchFunction(
    T *object, bool (T::*fkt)(Bead *, Bead *, Bead *, const Eigen::Vector3d &,
                              const Eigen::Vector3d &, const Eigen::Vector3d &,
                              const double dist12, const double dist13,
                              const double dist23)) {
  if (_match_function) {
    delete _match_function;
  }
  _match_function = dynamic_cast<Functor *>(new FunctorMember<T>(object, fkt));
}

inline void NBList_3Body::SetMatchFunction(bool (*fkt)(
    Bead *, Bead *, Bead *, const Eigen::Vector3d &, const Eigen::Vector3d &,
    const Eigen::Vector3d &, const double dist12, const double dist13,
    const double dist23)) {
  if (_match_function) {
    delete _match_function;
  }
  _match_function = dynamic_cast<Functor *>(new FunctorNonMember(fkt));
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_NBLIST_3BODY_H */
