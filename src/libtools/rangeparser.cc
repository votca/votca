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

#include <functional>
#include <iostream>
#include <stdexcept>
#include <votca/tools/rangeparser.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace tools {

RangeParser::RangeParser()
    //: _has_begin(false) , _has_end(false)
    = default;

void RangeParser::Parse(std::string str) {
  // remove all spaces in string
  std::string::iterator it = std::remove_if(
      str.begin(), str.end(), std::bind2nd(std::equal_to<char>(), ' '));
  str = std::string(str.begin(), it);

  Tokenizer tok(str, ",");
  for (std::string bl : tok) ParseBlock(bl);
}

void RangeParser::ParseBlock(std::string str) {
  Tokenizer tokenizer(str, ":");
  std::vector<std::string> toks;

  block_t block;
  block._stride = 1;

  tokenizer.ToVector(toks);
  if (toks.size() > 3 || toks.size() < 1) {
    throw std::runtime_error("invalid range");
  }

  block._begin = block._end = std::stoi(toks[0]);

  if (toks.size() == 2) block._end = std::stoi(toks[1]);

  if (toks.size() == 3) {
    block._stride = std::stoi(toks[1]);
    block._end = std::stoi(toks[2]);
  }

  if (block._begin * block._stride > block._end * block._stride) {
    throw std::runtime_error(
        std::string("invalid range " + str +
                    ": begin, end and stride do not form a closed interval"));
  }

  _blocks.push_back(block);
}

RangeParser::iterator& RangeParser::iterator::operator++() {
  _current += (*_block)._stride;
  if (_current > (*_block)._end) {
    ++_block;
    if (_block != _parent->_blocks.end())
      _current = (*_block)._begin;
    else
      _current = -1;
  }
  return *this;
}

}  // namespace tools
}  // namespace votca
