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

#ifndef VOTCA_TOOLS_RANGEPARSER_H
#define VOTCA_TOOLS_RANGEPARSER_H

// Standard includes
#include <list>
#include <ostream>
#include <string>

// Local VOTCA includes
#include "types.h"

namespace votca {
namespace tools {

/**
 * \brief RangeParser
 *
 * parse strings like min:step:max, not flexible enough yet to be really useful
 */
class RangeParser {
 public:
  RangeParser();

  void Parse(std::string str);

  void Add(Index begin, Index end, Index stride = 1);

 private:
  struct block_t {
    block_t() = default;
    block_t(const Index &begin, const Index &end, const Index &stride)
        : begin_(begin), end_(end), stride_(stride) {}

    Index begin_, end_, stride_;
  };

 public:
  struct iterator {
    iterator() = default;

    Index operator*() const { return current_; }

    RangeParser::iterator &operator++();

    bool operator==(const RangeParser::iterator &);
    bool operator!=(const RangeParser::iterator &);

   private:
    RangeParser *parent_;

    iterator(RangeParser *, std::list<block_t>::iterator);
    std::list<block_t>::iterator block_;
    Index current_;

    friend class RangeParser;
  };

  RangeParser::iterator begin();
  RangeParser::iterator end();

 private:
  void ParseBlock(std::string str);

  std::list<block_t> blocks_;

  friend std::ostream &operator<<(std::ostream &out, const RangeParser &rp);
};

inline void RangeParser::Add(Index begin, Index end, Index stride) {
  blocks_.push_back(block_t(begin, end, stride));
}

inline RangeParser::iterator RangeParser::begin() {
  return RangeParser::iterator(this, blocks_.begin());
}

inline RangeParser::iterator RangeParser::end() {
  return RangeParser::iterator(this, blocks_.end());
}

inline RangeParser::iterator::iterator(RangeParser *parent,
                                       std::list<block_t>::iterator block)
    : parent_(parent), block_(block) {
  if (block != parent->blocks_.end()) {
    current_ = (*block).begin_;
  } else {
    current_ = -1;
  }
}

inline bool RangeParser::iterator::operator==(const RangeParser::iterator &i) {
  return (block_ == i.block_) && (current_ == i.current_);
}

inline bool RangeParser::iterator::operator!=(const RangeParser::iterator &i) {
  return !((block_ == i.block_) && (current_ == i.current_));
}

inline std::ostream &operator<<(std::ostream &out, const RangeParser &rp) {
  std::list<RangeParser::block_t>::const_iterator iter(rp.blocks_.begin());
  for (; iter != rp.blocks_.end(); ++iter) {
    if (iter != rp.blocks_.begin()) {
      out << ",";
    }
    if (iter->begin_ == iter->end_) {
      out << iter->begin_;
    } else if (iter->stride_ == 1) {
      out << iter->begin_ << ":" << iter->end_;
    } else {
      out << iter->begin_ << ":" << iter->stride_ << ":" << iter->end_;
    }
  }
  return out;
}

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_RANGEPARSER_H
