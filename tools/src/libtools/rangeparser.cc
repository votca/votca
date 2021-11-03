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

// Standard includes
#include <functional>
#include <iostream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/tools/rangeparser.h"
#include "votca/tools/tokenizer.h"

namespace votca {
namespace tools {

RangeParser::RangeParser() = default;

void RangeParser::Parse(std::string str) {
  // remove all spaces in string
  std::string::iterator it =
      std::remove_if(str.begin(), str.end(), [](char a) { return a == ' '; });
  str = std::string(str.begin(), it);

  Tokenizer tok(str, ",");
  for (std::string bl : tok) {
    ParseBlock(bl);
  }
}

void RangeParser::ParseBlock(std::string str) {
  std::vector<std::string> toks = Tokenizer(str, ":").ToVector();

  block_t block;
  block.stride_ = 1;

  if (toks.size() > 3 || toks.size() < 1) {
    throw std::runtime_error("invalid range");
  }

  block.begin_ = block.end_ = std::stoi(toks[0]);

  if (toks.size() == 2) {
    block.end_ = std::stoi(toks[1]);
  }

  if (toks.size() == 3) {
    block.stride_ = std::stoi(toks[1]);
    block.end_ = std::stoi(toks[2]);
  }

  if (block.begin_ * block.stride_ > block.end_ * block.stride_) {
    throw std::runtime_error(
        std::string("invalid range " + str +
                    ": begin, end and stride do not form a closed interval"));
  }

  blocks_.push_back(block);
}

RangeParser::iterator& RangeParser::iterator::operator++() {
  current_ += (*block_).stride_;
  if (current_ > (*block_).end_) {
    ++block_;
    if (block_ != parent_->blocks_.end()) {
      current_ = (*block_).begin_;
    } else {
      current_ = -1;
    }
  }
  return *this;
}

}  // namespace tools
}  // namespace votca
